#!/usr/bin/env python3
"""BRAKER3 input handoff, command construction and prediction integrity checks."""
import argparse
import base64
import hashlib
import json
import re
import shlex
import subprocess
from pathlib import Path
from urllib.parse import unquote

from reference_repeats import fasta, sha
from prepare_protein_evidence import records


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def decode(value):
    return json.loads(base64.b64decode(value))


def plan(preflight, repeats, proteins, rna):
    """Select matched, ready evidence; never silently downgrade a prediction mode."""
    for label, data in (('preflight', preflight), ('reference/repeats', repeats),
                        ('proteins', proteins), ('RNA alignment', rna)):
        if data.get('status') != 'ready':
            raise ValueError(f'{label} is not ready')
    refs = {r['assembly_id']: r for r in repeats['assemblies']}
    prots = {r['assembly_id']: r for r in proteins['assemblies']}
    alignments = {(r['assembly_id'], r['readset_id']): r for r in rna['readsets']}
    requested_reads = {r['readset_id']: r for r in preflight['rna_readsets']}
    bam_paths = {}
    for library in rna['libraries']:
        if library['status'] != 'ready':
            raise ValueError('RNA library is not ready')
        if len(library['readsets']) != len(library['bams']) or len(library['bams']) != len(library['bam_indices']):
            raise ValueError('RNA library BAM/index count differs from readsets')
        for rid, bam, index in zip(library['readsets'], library['bams'], library['bam_indices']):
            key = (library['assembly_id'], rid)
            if key in bam_paths:
                raise ValueError('Duplicate RNA library readset')
            bam_paths[key] = (bam, index)
    jobs = []
    for row in preflight['assemblies']:
        aid = row['assembly_id']
        if int(row['settings']['genetic_code']) != 1:
            raise ValueError(f'{aid}: this BRAKER3 stage supports the standard genetic code only')
        if row['prediction_mode'] != 'rna_protein':
            raise ValueError(f'{aid}: this BRAKER3 stage requires RNA and proteins; mode {row["prediction_mode"]} is not implemented yet')
        ref, prot = refs[aid], prots[aid]
        if ref['status'] != 'ready' or prot['status'] != 'ready' or prot['protein_count'] < 1:
            raise ValueError(f'{aid}: reference/repeat or protein evidence is not ready')
        rid_list = sorted(row['rna_readsets'])
        if not rid_list or len(rid_list) != len(set(rid_list)):
            raise ValueError(f'{aid}: missing or duplicate RNA readsets')
        selected = []
        for rid in rid_list:
            alignment = alignments[(aid, rid)]
            current = requested_reads[rid]
            if any(current.get(key) != alignment.get(key) for key in ('library_id', 'sample_id', 'layout')):
                raise ValueError(f'{aid}/{rid}: current RNA sheet differs from accepted BAM metadata')
            strand = alignment['strand']
            if alignment['status'] != 'ready' or strand['status'] != 'ready':
                raise ValueError(f'{aid}/{rid}: BAM or strand inference is not ready')
            if strand['strandedness'] not in ('forward', 'reverse', 'unstranded'):
                raise ValueError(f'{aid}/{rid}: unresolved strandedness')
            if current['strandedness'] != 'auto' and current['strandedness'] != strand['strandedness']:
                raise ValueError(f'{aid}/{rid}: current strand setting differs from accepted RNA QC; rerun alignment reporting')
            if alignment['reference_identity'] != ref['repeats']['reference_identity']:
                raise ValueError(f'{aid}/{rid}: alignment and masked reference identities differ')
            if alignment.get('spliced_primary_reads', 0) < 1:
                raise ValueError(f'{aid}/{rid}: no spliced RNA evidence')
            selected.append({key: alignment[key] for key in ('readset_id', 'library_id', 'sample_id', 'layout')})
            selected[-1]['strandedness'] = strand['strandedness']
            selected[-1]['bam'], selected[-1]['bam_index'] = bam_paths[(aid, rid)]
        if len({r['bam'] for r in selected}) != len(selected):
            raise ValueError(f'{aid}: multiple readsets reference the same BAM')
        libraries = {}
        for readset in selected:
            prior = libraries.setdefault(readset['library_id'], readset['strandedness'])
            if prior != readset['strandedness']:
                raise ValueError(f'{aid}: conflicting lane strands within a library')
        token = hashlib.sha256(aid.encode()).hexdigest()[:10]
        jobs.append(dict(assembly_id=aid, mode='rna_protein',
                         species=re.sub(r'[^A-Za-z0-9_]', '_', aid) + '_braker_' + token,
                         genetic_code=row['settings']['genetic_code'],
                         reference_identity=ref['repeats']['reference_identity'],
                         masked_sha256=ref['repeats']['masked_sha256'],
                         proteins_sha256=prot['proteins_sha256'], protein_count=prot['protein_count'],
                         protein_source_ids=prot['source_ids'], readsets=selected))
    return {'schema_version': 1, 'status': 'ready', 'jobs': jobs}


def validate_inputs(spec, genome, proteins, bams, indices, identity):
    import pysam
    if not bams or len(bams) != len(spec['readsets']) or len(indices) != len(bams):
        raise ValueError('BRAKER BAM/index/readset count mismatch')
    if sha(genome) != spec['masked_sha256'] or sha(proteins) != spec['proteins_sha256']:
        raise ValueError('BRAKER genome/protein bytes differ from accepted evidence')
    expected = identity['contigs']
    if identity['reference_identity'] != spec['reference_identity']:
        raise ValueError('BRAKER reference identity differs from accepted BAM reference')
    observed = [{'id': h.split()[0], 'length': len(s),
                 'uppercase_sha256': hashlib.sha256(s.upper().encode()).hexdigest()}
                for h, s in fasta(genome)]
    if observed != expected:
        raise ValueError('BRAKER softmasked genome changed contig IDs/order/lengths/bases')
    dictionary = [(r['id'], r['length']) for r in expected]
    for path, index, readset in zip(bams, indices, spec['readsets']):
        if not Path(index).is_file():
            raise ValueError('BRAKER BAM index missing')
        pysam.quickcheck(path)
        with pysam.AlignmentFile(path, 'rb', index_filename=str(index)) as bam:
            if list(zip(bam.references, bam.lengths)) != dictionary or bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
                raise ValueError('BRAKER BAM dictionary/order differs from reference')
            bam.check_index()
            groups = {r['ID']: r for r in bam.header.to_dict().get('RG', [])}
            group = groups.get(readset['readset_id'], {})
            if group.get('LB') != readset['library_id'] or group.get('SM') != (readset.get('sample_id') or readset['library_id']):
                raise ValueError('BRAKER BAM is associated with the wrong readset/library')
    return {**spec, 'status': 'ready', 'validation_scope': 'Accepted genome/protein checksums, reference identity, BAM headers and readable indices; upstream streaming BAM QC retained'}


def command(spec, genome, proteins, bams, threads, rounds, workdir, augustus, ab_initio=False, check=False):
    if spec['mode'] != 'rna_protein' or not bams:
        raise ValueError('BRAKER3 requires RNA and protein evidence')
    if not 1 <= threads <= 48 or rounds < 1:
        raise ValueError('BRAKER requires 1..48 tool threads and positive training rounds')
    # --stranded is for strand-separated coverage BAMs/UTRs, not library protocols.
    args = ['braker.pl', f'--threads={threads}', f'--genome={Path(genome).resolve()}',
            f'--prot_seq={Path(proteins).resolve()}',
            '--bam=' + ','.join(str(Path(p).resolve()) for p in bams),
            f'--species={spec["species"]}', '--gff3',
            f'--AUGUSTUS_CONFIG_PATH={Path(augustus).resolve()}',
            f'--workingdir={Path(workdir).resolve()}', f'--rounds={rounds}',
            f'--translation_table={spec["genetic_code"]}', '--verbosity=3']
    if ab_initio:
        args.append('--AUGUSTUS_ab_initio')
    if check:
        args.append('--checkSoftware')
    return args


def run(spec, request, output):
    output = Path(output)
    output.mkdir(exist_ok=True)
    version = subprocess.run(['braker.pl', '--version'], check=True, capture_output=True, text=True)
    text = version.stdout + version.stderr
    if not re.search(r'\bversion\s+3\.', text):
        raise ValueError('The configured runtime is not BRAKER3: ' + text)
    (output / 'braker_version.txt').write_text(text, encoding='utf-8')
    args = command(spec, request['genome'], request['proteins'], request['bams'],
                   request['threads'], request['rounds'], request['workdir'],
                   request['augustus'], request['ab_initio'], request['check'])
    dump(output / 'command.json', args)
    dump(output / 'inputs.json', {**spec, 'execution_mode': 'software_check' if request['check'] else 'prediction',
                                 'container': request['container'], 'threads': request['threads'],
                                 'rounds': request['rounds'], 'augustus_ab_initio': request['ab_initio'],
                                 'strand_handling': 'Splice-derived XS tags; no UTR prediction or strand-separated BAM flag'})
    print(shlex.join(args), flush=True)
    # Stream diagnostics into native Nextflow output, including failed scratch tasks.
    with (output / 'console.log').open('w', encoding='utf-8') as log:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='', flush=True)
            log.write(line)
        code = process.wait()
    if code:
        raise SystemExit(code)


def validate_outputs(spec, genome_identity, gff, proteins, cds):
    lengths = {r['id']: r['length'] for r in genome_identity['contigs']}
    ids, parents, coding, genes, transcripts = {}, [], set(), set(), set()
    features = 0
    with Path(gff).open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Prediction GFF3 line {number}: expected 9 fields')
            seq, _, kind, start, end, _, strand, phase, attributes = fields
            if seq not in lengths or not 1 <= int(start) <= int(end) <= lengths[seq]:
                raise ValueError('Prediction feature outside accepted reference')
            if strand not in ('+', '-') or (kind == 'CDS' and phase not in ('0', '1', '2')):
                raise ValueError('Prediction has invalid strand or CDS phase')
            attrs = dict(item.split('=', 1) for item in attributes.split(';') if '=' in item)
            identifier = unquote(attrs['ID']) if 'ID' in attrs else None
            if identifier:
                if identifier in ids and not (kind == 'CDS' and ids[identifier] == 'CDS'):
                    raise ValueError('Duplicate prediction feature ID')
                ids[identifier] = kind
            if kind in ('gene', 'mRNA', 'transcript') and not identifier:
                raise ValueError('Prediction gene/transcript is missing ID')
            if kind == 'gene':
                genes.add(identifier)
            if kind in ('mRNA', 'transcript'):
                transcripts.add(identifier)
            relations = [unquote(x) for x in attrs.get('Parent', '').split(',') if x]
            if kind != 'gene' and not relations:
                raise ValueError('Prediction child feature is missing Parent')
            parents.extend(relations)
            if kind == 'CDS':
                coding.update(relations)
            features += 1
    if not genes or not coding or not set(parents) <= ids.keys() or not coding <= transcripts:
        raise ValueError('Prediction is empty or has unresolved feature parents')
    def sequences(path):
        result = {}
        with Path(path).open(encoding='utf-8') as handle:
            for header, sequence in records(handle):
                name = header.split()[0]
                if name in result or not sequence:
                    raise ValueError('Duplicate/empty predicted sequence')
                result[name] = len(sequence)
        return result
    aa, nucleotide = sequences(proteins), sequences(cds)
    if set(aa) != coding or set(nucleotide) != coding:
        raise ValueError('Predicted protein/CDS IDs differ from coding transcripts')
    return {'assembly_id': spec['assembly_id'], 'status': 'ready', 'genes': len(genes),
            'coding_transcripts': len(coding), 'features': features, 'proteins': len(aa),
            'cds_sequences': len(nucleotide), 'reference_identity': spec['reference_identity'],
            'gff3_sha256': sha(gff), 'proteins_sha256': sha(proteins), 'cds_sha256': sha(cds),
            'validation_scope': 'Nonempty predictions, reference bounds, feature IDs/parents and matching protein/CDS IDs; biological quality remains to be assessed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan', 'validate-inputs', 'run', 'validate-outputs'])
    for key in ('preflight', 'repeats', 'protein-status', 'rna-status', 'spec', 'request', 'output', 'identity', 'genome', 'proteins', 'gff', 'cds'):
        parser.add_argument('--' + key)
    args = parser.parse_args()
    if args.action == 'plan':
        dump(args.output, plan(load(args.preflight), load(args.repeats), load(args.protein_status), load(args.rna_status)))
    elif args.action == 'run':
        run(decode(args.spec), decode(args.request), args.output)
    elif args.action == 'validate-inputs':
        request = decode(args.request)
        dump(args.output, validate_inputs(decode(args.spec), request['genome'], request['proteins'],
                                         request['bams'], request['indices'], load(args.identity)))
    else:
        dump(args.output, validate_outputs(decode(args.spec), load(args.identity), args.gff, args.proteins, args.cds))


if __name__ == '__main__':
    main()
