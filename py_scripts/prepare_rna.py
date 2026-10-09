#!/usr/bin/env python3
"""Assembly-independent short-read acquisition, integrity checks, trimming and subsets."""
import argparse
import base64
import gzip
import hashlib
import itertools
import json
import math
import shutil
import subprocess
from pathlib import Path


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def compressed(path):
    with Path(path).open('rb') as handle:
        return handle.read(2) == b'\x1f\x8b'


def records(path):
    opener = gzip.open if compressed(path) else open
    with opener(path, 'rt', encoding='ascii') as handle:
        for number in itertools.count(1):
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = (handle.readline() for _ in range(3))
            lines = [x.rstrip('\r\n') for x in (header, sequence, plus, quality)]
            header, sequence, plus, quality = lines
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'{path}: invalid four-line FASTQ record {number}')
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f'{path}: sequence/quality length mismatch at record {number}')
            if not set(sequence.upper()) <= set('ACGTRYSWKMBDHVN'):
                raise ValueError(f'{path}: invalid nucleotide at record {number}')
            if any(not 33 <= ord(x) <= 126 for x in quality):
                raise ValueError(f'{path}: invalid quality encoding at record {number}')
            if len(header.split()[0]) <= 1:
                raise ValueError(f'{path}: empty read ID at record {number}')
            yield tuple(lines)


def read_id(record):
    token = record[0].split()[0][1:]
    return token[:-2] if token.endswith(('/1', '/2')) else token


def pairs(files):
    streams = [records(path) for path in files]
    for index, group in enumerate(itertools.zip_longest(*streams), 1):
        if any(record is None for record in group):
            raise ValueError(f'Mates have different read counts at pair {index}')
        if len(group) == 2:
            if read_id(group[0]) != read_id(group[1]):
                raise ValueError(f'Mate IDs disagree at pair {index}')
            for mate, record in enumerate(group, 1):
                fields = record[0].split()
                if fields[0].endswith(('/1', '/2')) and not fields[0].endswith(f'/{mate}'):
                    raise ValueError(f'Mate suffix disagrees with input order at pair {index}')
                if len(fields) > 1 and fields[1].startswith(('1:', '2:')) and not fields[1].startswith(f'{mate}:'):
                    raise ValueError(f'CASAVA mate field disagrees with input order at pair {index}')
        yield group


def inspect(files, layout):
    if len(files) != (2 if layout == 'paired' else 1):
        raise ValueError('Actual file count does not match declared layout')
    count, bases = 0, [0] * len(files)
    for group in pairs(files):
        count += 1
        for mate, record in enumerate(group):
            bases[mate] += len(record[1])
    if not count:
        raise ValueError('No usable reads remain')
    return {'read_units': count, 'read_count': count * len(files), 'bases_by_mate': bases,
            'files': [Path(x).name for x in files], 'sha256': [digest(x) for x in files]}


def plan(manifest, options, outdir):
    if manifest['status'] != 'ready':
        raise ValueError('RNA preparation requires ready input preflight')
    used = {rid for assembly in manifest['assemblies'] for rid in assembly['rna_readsets']}
    rows = sorted((row for row in manifest['rna_readsets'] if row['readset_id'] in used),
                  key=lambda row: row['readset_id'])
    quality = int(options['quality'])
    window, minimum, max_n = int(options['window']), int(options['min_length']), int(options['max_n'])
    if not 0 <= quality <= 93 or window < 1 or minimum < 1 or max_n < 0:
        raise ValueError('Invalid RNA trimming thresholds')
    for row in rows:
        row['raw_request'] = {key: row.get(key) for key in
                              ('readset_id', 'layout', 'fastq_1', 'fastq_2', 'run_accession')}
        row['trim_request'] = {key: row[key] for key in ('trim_adapters', 'trim_quality')}
        row['trim_request'].update(quality=quality, window=window, min_length=minimum, max_n=max_n)
        row['subset_request'] = {key: row[key] for key in ('subset_fraction', 'subset_seed')}
    result = {'schema_version': 1, 'readsets': rows,
              'assemblies': [{'assembly_id': row['assembly_id'], 'rna_readsets': row['rna_readsets']}
                             for row in manifest['assemblies']]}
    Path(outdir).mkdir()
    write_json(Path(outdir) / 'plan.json', result)


def tool(command, directory, log_name):
    with (directory / log_name).open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def bundle_result(spec, outdir, stage):
    Path(outdir).mkdir()
    return {'schema_version': 1, 'readset_id': spec['readset_id'], 'stage': stage,
            'status': 'failed', 'request': spec, 'implementation_sha256': digest(__file__)}


def acquire(spec, inputs, outdir, threads):
    outdir = Path(outdir)
    result = bundle_result(spec, outdir, 'raw')
    try:
        accession = spec.get('run_accession')
        if accession:
            result['sra_tools_version'] = subprocess.run(['fasterq-dump', '--version'], check=True,
                capture_output=True, text=True).stdout.strip()
            sra, extracted = outdir / 'sra', outdir / 'extracted'
            extracted.mkdir()
            command = ['prefetch', accession, '--max-size', 'u', '-O', str(sra)]
            write_json(outdir / 'prefetch_command.json', command)
            tool(command, outdir, 'prefetch.log')
            command = ['fasterq-dump', str(sra / accession / (accession + '.sra')),
                       '--split-files', '-e', str(threads), '-t', str(outdir / 'scratch'),
                       '-O', str(extracted)]
            write_json(outdir / 'fasterq_command.json', command)
            tool(command, outdir, 'fasterq.log')
            inputs = sorted(extracted.glob('*.fastq'))
            expected = [extracted / f'{accession}_{mate}.fastq' for mate in (1, 2)]
            if spec['layout'] == 'paired' and inputs != expected:
                raise ValueError('Accession output does not match requested paired layout')
            if spec['layout'] == 'single' and (len(inputs) != 1 or inputs[0].name not in
                    (accession + '.fastq', accession + '_1.fastq')):
                raise ValueError('Accession output does not match requested single layout')
        files = []
        for mate, original in enumerate(inputs, 1):
            original = Path(original)
            if accession:
                target = outdir / f'R{mate}.fastq.gz'
                with original.open('rb') as source, gzip.open(target, 'wb') as destination:
                    shutil.copyfileobj(source, destination, 1024 * 1024)
            else:
                target = original  # Nextflow already staged this declared input.
            files.append(target)
        result.update(inspect(files, spec['layout']), layout=spec['layout'], status='ready',
                      quality_encoding='Phred+33 expected; printable-range validation only')
        if accession:
            # The uncompressed conversion and conversion scratch are transient.
            # Keep downloaded SRA and compressed raw reads for reuse/provenance.
            shutil.rmtree(extracted)
            scratch = outdir / 'scratch'
            if scratch.exists():
                shutil.rmtree(scratch)
    except (ValueError, OSError, EOFError, UnicodeError, subprocess.CalledProcessError) as exc:
        result['reason'] = str(exc)
    write_json(outdir / 'rna.json', result)
    return result


def metadata(path):
    path = Path(path)
    return read_json(path / 'rna.json' if path.is_dir() else path)


def fastp_command(raw, settings, outdir, threads, reads):
    data = metadata(raw)
    files = [Path(path) for path in reads]
    outputs = [Path(outdir) / f'R{mate}.fastq.gz' for mate in range(1, len(files) + 1)]
    command = ['fastp', '--in1', str(files[0]), '--out1', str(outputs[0]), '--thread', str(threads),
               '--json', str(Path(outdir) / (data['readset_id'] + '_fastp.json')),
               '--html', str(Path(outdir) / (data['readset_id'] + '_fastp.html')),
               '--report_title', data['readset_id'], '--disable_trim_poly_g',
               '--disable_quality_filtering', '--length_required', str(settings['min_length']),
               '--n_base_limit', str(settings['max_n'])]
    if len(files) == 2:
        command += ['--in2', str(files[1]), '--out2', str(outputs[1])]
        if settings['trim_adapters']:
            command.append('--detect_adapter_for_pe')
    if not settings['trim_adapters']:
        command.append('--disable_adapter_trimming')
    if settings['trim_quality']:
        command += ['--cut_tail', '--cut_tail_window_size', str(settings['window']),
                    '--cut_tail_mean_quality', str(settings['quality'])]
    return command, outputs


def trim(rid, settings, raw, outdir, threads, reads):
    outdir = Path(outdir)
    result = bundle_result({'readset_id': rid, **settings}, outdir, 'processed')
    try:
        data = metadata(raw)
        if not settings['trim_adapters'] and not settings['trim_quality']:
            result.update({key: data[key] for key in
                           ('read_units', 'read_count', 'bases_by_mate', 'sha256', 'layout')},
                          files=[Path(path).name for path in reads], passthrough=True,
                          status='ready', raw_read_units=data['read_units'], raw_sha256=data['sha256'],
                          retained_fraction=1.0, processing='passthrough: trimming disabled')
            write_json(outdir / 'rna.json', result)
            return result
        command, files = fastp_command(raw, settings, outdir, threads, reads)
        write_json(outdir / 'command.json', command)
        version = subprocess.run(['fastp', '--version'], check=True, capture_output=True, text=True)
        result['fastp_version'] = (version.stdout + version.stderr).strip()
        tool(command, outdir, 'fastp.log')
        result.update(inspect(files, data['layout']), layout=data['layout'], status='ready',
                      raw_read_units=data['read_units'], raw_sha256=data['sha256'])
        result['retained_fraction'] = result['read_units'] / data['read_units']
    except (ValueError, OSError, EOFError, UnicodeError, subprocess.CalledProcessError) as exc:
        result['reason'] = str(exc)
    write_json(outdir / 'rna.json', result)
    return result


def subset(rid, settings, clean, outdir, reads):
    outdir = Path(outdir)
    result = bundle_result({'readset_id': rid, **settings}, outdir, 'selected')
    try:
        data = metadata(clean)
        inputs = [Path(path) for path in reads]
        fraction, seed = float(settings['subset_fraction']), int(settings['subset_seed'])
        if not math.isfinite(fraction) or not 0 < fraction <= 1 or seed < 0:
            raise ValueError('Invalid RNA subset settings')
        if fraction == 1:
            result['passthrough'] = True
            stats = {key: data[key] for key in ('read_units', 'read_count', 'bases_by_mate', 'sha256')}
            stats['files'] = [path.name for path in inputs]
        else:
            outputs = [outdir / f'R{mate}.fastq.gz' for mate in range(1, len(inputs) + 1)]
            # Select both mates together by read ID; stable across file ordering and runs.
            handles = [gzip.open(path, 'wt', encoding='ascii') for path in outputs]
            try:
                threshold = int(fraction * (1 << 64))
                for group in pairs(inputs):
                    value = int.from_bytes(hashlib.sha256(f'{seed}\0{read_id(group[0])}'.encode()).digest()[:8], 'big')
                    if value < threshold:
                        for handle, record in zip(handles, group):
                            handle.write('\n'.join(record) + '\n')
            finally:
                for handle in handles:
                    handle.close()
            stats = inspect(outputs, data['layout'])
        result.update(stats, layout=data['layout'], status='ready', clean_read_units=data['read_units'],
                      clean_sha256=data['sha256'], selection='sha256(seed + NUL + read_id), first 64 bits')
        result['processing'] = {key: data.get(key) for key in
                                ('request', 'raw_read_units', 'read_units', 'retained_fraction', 'fastp_version')}
    except (ValueError, OSError, EOFError, UnicodeError) as exc:
        result['reason'] = str(exc)
    write_json(outdir / 'rna.json', result)
    return result


def qc(rid, stage, reads, outdir, threads):
    outdir = Path(outdir)
    outdir.mkdir()
    # Nextflow supplies readset/stage-prefixed input names through stageAs.
    command = ['fastqc', '--threads', str(threads), '--outdir', str(outdir), *map(str, reads)]
    write_json(outdir / 'command.json', command)
    version = subprocess.run(['fastqc', '--version'], check=True, capture_output=True, text=True)
    (outdir / 'version.txt').write_text(version.stdout + version.stderr)
    tool(command, outdir, 'fastqc.log')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan', 'acquire', 'trim', 'subset', 'qc'])
    parser.add_argument('--request', required=True)
    parser.add_argument('--manifest')
    parser.add_argument('--inputs', default='W10=')
    parser.add_argument('--bundle')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    spec = json.loads(base64.b64decode(args.request))
    inputs = json.loads(base64.b64decode(args.inputs))
    if args.action == 'plan':
        plan(read_json(args.manifest), spec, args.outdir)
    elif args.action == 'acquire':
        acquire(spec, inputs, args.outdir, args.threads)
    elif args.action == 'trim':
        trim(spec['readset_id'], spec, args.bundle, args.outdir, args.threads, inputs)
    elif args.action == 'subset':
        subset(spec['readset_id'], spec, args.bundle, args.outdir, inputs)
    else:
        qc(spec['readset_id'], spec['stage'], inputs, args.outdir, args.threads)


if __name__ == '__main__':
    main()
