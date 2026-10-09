#!/usr/bin/env python3
"""Plan, acquire, validate and combine source-preserving prediction proteins."""
import argparse
import base64
import csv
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = 1
AA = set('ACDEFGHIKLMNPQRSTVWYBXZJUO')


def checksum(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def plan(manifest, catalogue, options, outdir):
    if manifest['status'] != 'ready':
        raise ValueError('Protein preparation requires a ready input preflight')
    outdir = Path(outdir)
    (outdir / 'sources').mkdir(parents=True)
    (outdir / 'assemblies').mkdir()
    sources, assemblies, local_checksums = {}, [], {}

    def local_checksum(path):
        key = str(Path(path).resolve())
        if key not in local_checksums:
            local_checksums[key] = checksum(key)
        return local_checksums[key]

    def register(spec):
        identity = {k: v for k, v in spec.items() if k not in ('local_path', 'local_sha256', 'retry_nonce')}
        if spec['kind'] == 'local':
            identity.update(local_path=spec['local_path'], local_sha256=spec['local_sha256'])
        digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        spec['source_id'] = 'p' + digest[:24]
        if spec['source_id'] in sources:
            previous = {k: v for k, v in sources[spec['source_id']].items() if k != 'retry_nonce'}
            if previous != spec:
                raise ValueError('Inconsistent source identity')
            return spec['source_id']
        published = options.get('published_sources')
        if published:
            previous_status = Path(published) / spec['source_id'] / 'source.json'
            if previous_status.is_file() and read_json(previous_status)['status'] == 'failed':
                # Failed tasks emit status data. Give only these sources a fresh
                # task identity on resume; successfully acquired sources stay cached.
                spec['retry_nonce'] = datetime.now(timezone.utc).isoformat()
        sources[spec['source_id']] = spec
        return spec['source_id']

    for row in manifest['assemblies']:
        settings = row['settings']
        ids = []
        if settings['use_protein_evidence']:
            clade = settings['orthodb_clade']
            if clade != 'none':
                if clade not in catalogue['partitions']:
                    raise ValueError(f'No pinned OrthoDB archive checksum for {clade}')
                spec = {'schema_version': SCHEMA, 'kind': 'orthodb', 'release': catalogue['release'],
                        'clade': clade, 'url': catalogue['catalogue_url'] + clade + '.fa.gz',
                        'expected_md5': catalogue['partitions'][clade]}
                cache = options.get('orthodb_database')
                if cache:
                    base = Path(cache).expanduser()
                    if not base.is_absolute():
                        base = Path(options['launch_dir']) / base
                    cached = base / (clade + '.fa.gz')
                    if cached.is_file():
                        spec.update(local_path=str(cached.resolve()), local_sha256=local_checksum(cached))
                ids.append(register(spec))
            for accession in settings['ref_proteome_accession']:
                ids.append(register({'schema_version': SCHEMA, 'kind': 'ncbi',
                                     'accession': accession, 'datasets_version': '18.35.0'}))
            for local in settings['protein_fastas']:
                ids.append(register({'schema_version': SCHEMA, 'kind': 'local',
                    'local_path': local, 'local_sha256': local_checksum(local)}))
            if not ids:
                raise ValueError(f"{row['assembly_id']}: enabled proteins have no sources")
        assembly = {'schema_version': SCHEMA, 'assembly_id': row['assembly_id'],
                    'enabled': settings['use_protein_evidence'], 'source_ids': sorted(set(ids))}
        assemblies.append(assembly)
        write_json(outdir / 'assemblies' / (row['assembly_id'] + '.json'), assembly)
    for sid, spec in sorted(sources.items()):
        write_json(outdir / 'sources' / (sid + '.json'), spec)
    result = {'schema_version': SCHEMA, 'sources': sorted(sources.values(), key=lambda x: x['source_id']),
              'assemblies': assemblies}
    write_json(outdir / 'plan.json', result)
    return result


def records(handle):
    header, sequence = None, []
    for line_number, line in enumerate(handle, 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence)
            header, sequence = line[1:].strip(), []
            if not header:
                raise ValueError(f'Empty FASTA header on line {line_number}')
        else:
            if header is None:
                raise ValueError(f'Sequence before first FASTA header on line {line_number}')
            sequence.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(sequence)


def normalize(handle, source_id, outdir):
    seen, count, residues, terminal_stops = set(), 0, 0, 0
    protein_path, crosswalk_path = outdir / 'proteins.fa', outdir / 'id_map.tsv'
    try:
        with protein_path.open('w', encoding='utf-8') as fasta, crosswalk_path.open('w', newline='', encoding='utf-8') as crosswalk:
            writer = csv.writer(crosswalk, delimiter='\t')
            writer.writerow(['protein_id', 'source_id', 'original_id', 'original_header',
                             'length_aa', 'terminal_stop_removed'])
            for header, sequence in records(handle):
                original_id = header.split()[0]
                if original_id in seen:
                    raise ValueError(f'Duplicate protein ID within source: {original_id}')
                seen.add(original_id)
                sequence = sequence.upper()
                removed_stop = sequence.endswith('*')
                if removed_stop:
                    sequence = sequence[:-1]
                if not sequence:
                    raise ValueError(f'Empty protein sequence: {original_id}')
                invalid = set(sequence) - AA
                if invalid:
                    raise ValueError(f'Invalid amino-acid characters/internal stop in {original_id}: {sorted(invalid)}')
                count += 1
                pid = f'{source_id}__{count:09d}'
                fasta.write('>' + pid + '\n')
                for offset in range(0, len(sequence), 60):
                    fasta.write(sequence[offset:offset + 60] + '\n')
                writer.writerow([pid, source_id, original_id, header, len(sequence), int(removed_stop)])
                residues += len(sequence)
                terminal_stops += int(removed_stop)
        if not count:
            raise ValueError('Protein FASTA contains no records')
    except Exception:
        # Never leave a partial FASTA looking like valid evidence.
        protein_path.unlink(missing_ok=True)
        crosswalk_path.unlink(missing_ok=True)
        raise
    return {'protein_count': count, 'residues': residues, 'terminal_stops_removed': terminal_stops,
            'proteins_sha256': checksum(protein_path), 'id_map_sha256': checksum(crosswalk_path)}


def download(url, destination):
    partial = destination.with_suffix(destination.suffix + '.partial')
    with urllib.request.urlopen(url, timeout=120) as response, partial.open('wb') as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.rename(destination)


def acquire(spec, staged_inputs, outdir):
    outdir = Path(outdir)
    outdir.mkdir()
    write_json(outdir / 'request.json', spec)
    result = {'schema_version': SCHEMA, 'source_id': spec['source_id'], 'request': spec,
              'status': 'failed', 'python_version': sys.version.split()[0],
              'implementation_sha256': checksum(__file__),
              'normalization': 'source namespace + record number; uppercase; remove one terminal stop; wrap 60',
              'acquired_at_utc': datetime.now(timezone.utc).isoformat()}
    try:
        local = None
        if spec.get('local_path'):
            if len(staged_inputs) != 1:
                raise ValueError('Expected exactly one staged source input')
            local = Path(staged_inputs[0])
            if checksum(local) != spec['local_sha256']:
                raise ValueError('Local source changed after planning; rerun input preparation')
        if spec['kind'] == 'orthodb':
            raw = outdir / 'raw.fa.gz'
            if local:
                shutil.copyfile(local, raw)
            else:
                download(spec['url'], raw)
            result.update(raw_file=raw.name, raw_sha256=checksum(raw),
                          archive_md5=checksum(raw, 'md5'), provider_md5=spec['expected_md5'])
            if result['archive_md5'] != spec['expected_md5']:
                raise ValueError('OrthoDB archive MD5 differs from the pinned provider checksum')
            with gzip.open(raw, 'rt', encoding='utf-8') as handle:
                result.update(normalize(handle, spec['source_id'], outdir))
        elif spec['kind'] == 'ncbi':
            version = subprocess.run(['datasets', 'version'], check=True, capture_output=True, text=True, timeout=60)
            result['datasets_version'] = version.stdout.strip()
            raw = outdir / 'raw.zip'
            command = ['datasets', 'download', 'genome', 'accession', spec['accession'],
                       '--include', 'protein', '--filename', str(raw), '--no-progressbar']
            write_json(outdir / 'download_command.json', command)
            with (outdir / 'download.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
            result.update(raw_file=raw.name, raw_sha256=checksum(raw))
            member = f"ncbi_dataset/data/{spec['accession']}/protein.faa"
            with zipfile.ZipFile(raw) as archive:
                proteins = [x for x in archive.namelist() if x.endswith('/protein.faa')]
                if proteins != [member]:
                    raise ValueError(f'Expected only exact-accession protein FASTA {member}; found {proteins}')
                with archive.open(member) as binary, io.TextIOWrapper(binary, encoding='utf-8') as handle:
                    result.update(normalize(handle, spec['source_id'], outdir))
                # Keep the archive (including its annotation metadata), not just a FASTA.
            result.update(protein_member=member)
        elif spec['kind'] == 'local':
            if local is None:
                raise ValueError('Missing staged local protein FASTA')
            with local.open('rb') as handle:
                compressed = handle.read(2) == b'\x1f\x8b'
            raw = outdir / ('raw.fa.gz' if compressed else 'raw.fa')
            shutil.copyfile(local, raw)
            result.update(raw_file=raw.name, raw_sha256=checksum(raw))
            opener = gzip.open if compressed else open
            with opener(raw, 'rt', encoding='utf-8') as handle:
                result.update(normalize(handle, spec['source_id'], outdir))
        else:
            raise ValueError(f"Unsupported source kind: {spec['kind']}")
        result['status'] = 'ready'
    except (ValueError, OSError, EOFError, UnicodeError, zipfile.BadZipFile,
            subprocess.CalledProcessError, subprocess.TimeoutExpired, urllib.error.URLError) as exc:
        result['reason'] = str(exc)
        (outdir / 'error.txt').write_text(str(exc) + '\n')
    write_json(outdir / 'source.json', result)
    return result


def combine(assembly, source_directory, outdir):
    outdir = Path(outdir)
    outdir.mkdir()
    result = {**assembly, 'status': 'not_requested', 'sources': [], 'protein_count': 0,
              'implementation_sha256': checksum(__file__)}
    if not assembly['enabled']:
        if assembly['source_ids']:
            raise ValueError('Disabled protein evidence has requested sources')
    else:
        bundles = {}
        for path in Path(source_directory).glob('*/source.json'):
            source = read_json(path)
            if source['source_id'] in bundles:
                raise ValueError('Duplicate source bundle')
            bundles[source['source_id']] = (path.parent, source)
        requested = assembly['source_ids']
        if not requested:
            raise ValueError('Enabled protein evidence has no sources')
        blocked = []
        for sid in requested:
            if sid not in bundles:
                blocked.append(f'{sid}: source bundle missing')
                continue
            folder, source = bundles[sid]
            result['sources'].append({'source_id': sid, 'status': source['status'],
                'request': source['request'], 'raw_sha256': source.get('raw_sha256'),
                'proteins_sha256': source.get('proteins_sha256'), 'protein_count': source.get('protein_count'),
                'reason': source.get('reason')})
            if source['status'] != 'ready':
                blocked.append(f"{sid}: {source.get('reason', 'source not ready')}")
            elif checksum(folder / 'proteins.fa') != source['proteins_sha256'] or checksum(folder / 'id_map.tsv') != source['id_map_sha256']:
                blocked.append(f'{sid}: source output checksum mismatch')
        if blocked:
            result.update(status='blocked', reason='; '.join(blocked))
        else:
            with (outdir / 'proteins.fa').open('wb') as proteins, (outdir / 'id_map.tsv').open('w', newline='', encoding='utf-8') as mapping:
                writer = csv.writer(mapping, delimiter='\t')
                for index, sid in enumerate(requested):
                    folder, source = bundles[sid]
                    with (folder / 'proteins.fa').open('rb') as handle:
                        shutil.copyfileobj(handle, proteins, length=1024 * 1024)
                    with (folder / 'id_map.tsv').open(newline='', encoding='utf-8') as handle:
                        reader = csv.reader(handle, delimiter='\t')
                        header = next(reader)
                        if index == 0:
                            writer.writerow(header)
                        writer.writerows(reader)
                    result['protein_count'] += source['protein_count']
            result.update(status='ready', proteins_sha256=checksum(outdir / 'proteins.fa'),
                          id_map_sha256=checksum(outdir / 'id_map.tsv'))
    write_json(outdir / 'evidence.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    planning = sub.add_parser('plan')
    planning.add_argument('--manifest', required=True)
    planning.add_argument('--catalogue', required=True)
    planning.add_argument('--options-base64', required=True)
    planning.add_argument('--outdir', required=True)
    sourcing = sub.add_parser('acquire')
    sourcing.add_argument('--spec', required=True)
    sourcing.add_argument('--staged-base64', required=True)
    sourcing.add_argument('--outdir', required=True)
    combining = sub.add_parser('combine')
    combining.add_argument('--assembly', required=True)
    combining.add_argument('--sources-directory', required=True)
    combining.add_argument('--outdir', required=True)
    args = parser.parse_args()
    decode = lambda value: json.loads(base64.b64decode(value).decode())
    if args.action == 'plan':
        plan(read_json(args.manifest), read_json(args.catalogue), decode(args.options_base64), args.outdir)
    elif args.action == 'acquire':
        acquire(read_json(args.spec), decode(args.staged_base64), args.outdir)
    elif args.action == 'combine':
        combine(read_json(args.assembly), args.sources_directory, args.outdir)


if __name__ == '__main__':
    main()
