#!/usr/bin/env python3
"""RNA alignment contracts and reports. Nextflow stages and transports every file."""
import argparse
import base64
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def decode(value):
    return json.loads(base64.b64decode(value))


def selected_reads(spec, paths):
    from prepare_rna import digest
    expected = spec.get('sha256') or []
    if len(paths) != (2 if spec['layout'] == 'paired' else 1) or len(expected) != len(paths):
        raise ValueError('Selected-read checksum/file count mismatch')
    if [digest(path) for path in paths] != expected:
        raise ValueError('FASTQ bytes differ from the accepted selected-read manifest; rerun RNA preparation')


def plan(manifest, prepared, locations):
    if manifest['status'] != 'ready' or prepared['status'] != 'ready':
        raise ValueError('Input preflight and selected RNA preparation must both be ready')
    selected = {row['readset_id']: row for row in prepared['readsets']}
    requested = {row['readset_id']: row for row in manifest['rna_readsets']}
    root = Path(locations['reference_dir'])
    references, jobs, assemblies = [], [], []
    for assembly in manifest['assemblies']:
        aid = assembly['assembly_id']
        assemblies.append({'assembly_id': aid, 'rna_readsets': assembly['rna_readsets']})
        if not assembly['rna_readsets']:
            continue
        reference = root / aid / (aid + '.filtered.fasta')
        identity = root / aid / (aid + '.contigs.json')
        data = load(identity)
        if data['assembly_id'] != aid or not reference.is_file():
            raise ValueError(f'Missing/mismatched normalized reference: {aid}')
        refid = 'ref_' + data['reference_identity']
        references.append({'assembly_id': aid, 'reference_id': refid,
                           'reference': str(reference.resolve()), 'identity': str(identity.resolve())})
        libraries = {}
        for rid in assembly['rna_readsets']:
            row, current = selected.get(rid), requested[rid]
            if not row or row['status'] != 'ready' or not row.get('read_paths'):
                raise ValueError(f'No ready prepared reads for {rid}')
            for key in ('library_id', 'sample_id', 'layout'):
                if row.get(key) != current.get(key):
                    raise ValueError(f'Prepared/current RNA identity differs: {rid}/{key}')
            signature = (row.get('sample_id'), row['layout'])
            lid = row['library_id']
            if lid in libraries and libraries[lid] != signature:
                raise ValueError(f'Conflicting sample/layout within library {aid}/{lid}')
            libraries[lid] = signature
            paths = [str(Path(x).resolve()) for x in row['read_paths']]
            if len(paths) != (2 if row['layout'] == 'paired' else 1) or not all(Path(x).is_file() for x in paths):
                raise ValueError(f'Prepared FASTQs missing for {rid}; retain their Nextflow work directory')
            spec = {key: row.get(key) for key in ('library_id', 'sample_id', 'layout', 'read_units', 'sha256')}
            spec.update(assembly_id=aid, readset_id=rid, reference_identity=data['reference_identity'])
            jobs.append({'assembly_id': aid, 'readset_id': rid, 'library_id': lid,
                         'reference_id': refid, 'read_paths': paths, 'spec': spec,
                         'policy': {'strandedness': current['strandedness'], **locations['strand_policy']}})
    return {'schema_version': 1, 'assemblies': assemblies, 'references': references, 'jobs': jobs}


def reference(source, identity):
    # Reuse the accepted reference parser and content identity contract.
    from reference_repeats import fasta, sha
    rows = [{'id': header.split()[0], 'length': len(seq),
             'uppercase_sha256': hashlib.sha256(seq.upper().encode('ascii')).hexdigest()}
            for header, seq in fasta(source)]
    if rows != identity['contigs']:
        raise ValueError('Indexed reference differs from the accepted contig identities')
    value = hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if value != identity['reference_identity']:
        raise ValueError('Invalid reference identity checksum')
    return {'contigs': rows, 'reference_identity': value, 'indexed_fasta_sha256': sha(source)}


def strand(qc, policy):
    result = {key: qc[key] for key in ('assembly_id', 'readset_id', 'library_id')}
    result.update(status='blocked', requested=policy['strandedness'], policy=policy,
                  informative=qc['strand_forward'] + qc['strand_reverse'])
    if qc['status'] != 'ready':
        result['reason'] = qc.get('reason', 'Alignment QC failed')
        return result
    requested = policy['strandedness']
    if requested not in ('auto', 'forward', 'reverse', 'unstranded'):
        raise ValueError('Invalid strandedness')
    if requested != 'auto':
        result.update(status='ready', strandedness=requested, method='explicit')
        return result
    n = result['informative']
    fraction = qc['strand_forward'] / n if n else None
    result['forward_fraction'] = fraction
    minimum, dominant, tolerance = int(policy['min_informative']), float(policy['dominant_fraction']), float(policy['unstranded_tolerance'])
    if minimum < 1 or not .5 < dominant <= 1 or not 0 <= tolerance < dominant - .5:
        raise ValueError('Invalid strand inference thresholds')
    if n < minimum:
        result['reason'] = f'Only {n} informative spliced reads; need {minimum}. Set strandedness explicitly if known.'
    elif fraction >= dominant:
        result.update(status='ready', strandedness='forward', method='splice_XS')
    elif qc['strand_reverse'] / n >= dominant:
        result.update(status='ready', strandedness='reverse', method='splice_XS')
    elif abs(fraction - .5) <= tolerance:
        result.update(status='ready', strandedness='unstranded', method='splice_XS')
    else:
        result['reason'] = 'Intermediate strand mixture; inspect libraries or provide a known strand setting'
    return result


def stats(bam_path, identity, spec, outdir):
    import pysam
    result = {**spec, 'status': 'blocked', 'strand_forward': 0, 'strand_reverse': 0}
    counts, junctions = Counter(), Counter()
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            expected = [(row['id'], row['length']) for row in identity['contigs']]
            if list(zip(bam.references, bam.lengths)) != expected:
                raise ValueError('BAM reference dictionary differs from indexed reference')
            if identity['reference_identity'] != spec['reference_identity']:
                raise ValueError('BAM/request reference identity mismatch')
            groups = {row['ID']: row for row in bam.header.to_dict().get('RG', [])}
            group = groups.get(spec['readset_id'], {})
            if group.get('LB') != spec['library_id'] or group.get('SM') != (spec.get('sample_id') or spec['library_id']):
                raise ValueError('BAM read group sample/library identity differs from request')
            previous, unmapped_tail = (-1, -1), False
            for read in bam.fetch(until_eof=True):
                if read.reference_id >= 0:
                    coordinate = (read.reference_id, read.reference_start)
                    if unmapped_tail or coordinate < previous:
                        raise ValueError('BAM is not coordinate sorted')
                    previous = coordinate
                else:
                    unmapped_tail = True
                if read.is_secondary or read.is_supplementary:
                    counts['secondary_or_supplementary'] += 1
                    continue
                counts['primary_reads'] += 1
                if not read.has_tag('RG') or read.get_tag('RG') != spec['readset_id']:
                    raise ValueError('Primary alignment has missing/incorrect read group')
                if (spec['layout'] == 'paired') != read.is_paired:
                    raise ValueError('BAM paired flag differs from requested layout')
                if read.is_paired:
                    if read.is_read1 == read.is_read2:
                        raise ValueError('Invalid primary mate flags')
                    counts['mate1_reads' if read.is_read1 else 'mate2_reads'] += 1
                if read.is_unmapped:
                    continue
                if read.reference_start < 0 or read.reference_end > bam.lengths[read.reference_id]:
                    raise ValueError('Alignment extends beyond reference')
                counts['mapped_reads'] += 1
                counts['properly_paired_reads'] += int(read.is_proper_pair)
                unique = read.has_tag('NH') and read.get_tag('NH') == 1
                counts['unique_mapped_reads'] += int(unique)
                counts['soft_clipped_reads'] += int(any(op == 4 for op, length in read.cigartuples or []))
                position, spliced = read.reference_start, False
                for op, length in read.cigartuples or []:
                    if op == 3:
                        spliced = True
                        if unique and read.mapping_quality >= spec['min_mapq']:
                            junctions[(read.reference_name, position + 1, position + length,
                                       read.get_tag('XS') if read.has_tag('XS') else '.')] += 1
                    if op in (0, 2, 3, 7, 8):
                        position += length
                counts['spliced_primary_reads'] += int(spliced)
                # Use one mate per fragment, avoiding double-counting pair orientations.
                if unique and spliced and read.mapping_quality >= spec['min_mapq'] and (not read.is_paired or read.is_read1) and read.has_tag('XS'):
                    xs = read.get_tag('XS')
                    if xs in ('+', '-'):
                        counts['strand_forward' if read.is_reverse == (xs == '-') else 'strand_reverse'] += 1
            expected_reads = int(spec['read_units']) * (2 if spec['layout'] == 'paired' else 1)
            if counts['primary_reads'] != expected_reads:
                raise ValueError(f'Expected {expected_reads} primary reads; found {counts["primary_reads"]}')
            if spec['layout'] == 'paired' and (counts['mate1_reads'] != spec['read_units'] or counts['mate2_reads'] != spec['read_units']):
                raise ValueError('Primary mate counts differ from selected read pairs')
            if not counts['mapped_reads']:
                raise ValueError('No mapped reads')
            result['status'] = 'ready'
    except (ValueError, OSError) as exc:
        result['reason'] = str(exc)
    result.update(counts)
    result['mapped_fraction'] = counts['mapped_reads'] / counts['primary_reads'] if counts['primary_reads'] else 0
    result['junctions'] = len(junctions)
    folder = Path(outdir)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / 'junctions.tsv').open('w', encoding='utf-8') as output:
        output.write('contig\tintron_start\tintron_end\tstrand\tsupporting_alignments\n')
        for key, count in sorted(junctions.items()):
            output.write('\t'.join(map(str, (*key, count))) + '\n')
    dump(folder / 'alignment.json', result)


def report(plan_data, statistics, strands, bams):
    qc = {(row['assembly_id'], row['readset_id']): row for row in statistics}
    resolved = {(row['assembly_id'], row['readset_id']): row for row in strands}
    aligned = {(row['assembly_id'], row['readset_id']): row for row in bams}
    libraries = {}
    for job in plan_data['jobs']:
        libraries.setdefault((job['assembly_id'], job['library_id']), []).append(job['readset_id'])
    rows = []
    for (aid, lid), ids in sorted(libraries.items()):
        states = [resolved.get((aid, rid), {}) for rid in ids]
        values = {row.get('strandedness') for row in states}
        valid = all(row.get('status') == 'ready' for row in states) and len(values) == 1
        paths = [aligned.get((aid, rid), {}) for rid in ids]
        row = dict(assembly_id=aid, library_id=lid, readsets=ids, status='blocked')
        if valid and all(paths):
            row.update(status='ready', strandedness=next(iter(values)), readsets=ids,
                       bams=[data['bam'] for data in paths], bam_indices=[data['bam_index'] for data in paths])
        else:
            row['reason'] = 'Missing readset BAM' if valid else 'Unresolved/conflicting readset strandedness'
        rows.append(row)
    assemblies = []
    for row in plan_data['assemblies']:
        relevant = [lib for lib in rows if lib['assembly_id'] == row['assembly_id']]
        assemblies.append({**row, 'status': ('not_requested' if not row['rna_readsets'] else
                           'ready' if relevant and all(lib['status'] == 'ready' for lib in relevant) else 'blocked'),
                           'bams': [path for lib in relevant if lib['status'] == 'ready' for path in lib['bams']]})
    return {'schema_version': 1, 'status': 'ready' if all(row['status'] != 'blocked' for row in assemblies) else 'incomplete',
            'assemblies': assemblies, 'libraries': rows,
            'readsets': [{**qc.get((job['assembly_id'], job['readset_id']), {}),
                          'strand': resolved.get((job['assembly_id'], job['readset_id']), {})} for job in plan_data['jobs']]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan', 'reads', 'reference', 'stats', 'strand', 'report'])
    for name in ('preflight', 'prepared', 'request', 'output', 'reference', 'identity', 'bam', 'outdir', 'plan', 'qc', 'strands'):
        parser.add_argument('--' + name)
    args = parser.parse_args()
    if args.action == 'reads':
        selected_reads(decode(args.request)['spec'], decode(args.request)['paths'])
    elif args.action == 'plan':
        dump(args.output, plan(load(args.preflight), load(args.prepared), decode(args.request)))
    elif args.action == 'reference':
        dump(args.output, reference(args.reference, load(args.identity)))
    elif args.action == 'stats':
        stats(args.bam, load(args.identity), decode(args.request), args.outdir)
    elif args.action == 'strand':
        dump(args.output, strand(load(args.identity), decode(args.request)))
    else:
        result = report(load(args.plan), [load(x) for x in Path(args.qc).glob('*.json')],
                        [load(x) for x in Path(args.strands).glob('*.json')],
                        decode(args.request))
        dump(Path(args.outdir) / 'rna_alignment_status.json', result)
        lines = ['# RNA alignment', '', f'Status: **{result["status"]}**', '',
                 'Mapping uses the accepted normalized reference and previously selected reads.',
                 'Automatic strandedness uses unique, high-quality spliced primary alignments (R1 for paired reads).',
                 'Validated readset BAMs are listed for BRAKER; gene prediction has not run.', '',
                 '| Assembly | Readset | Mapped | Informative | Strand | Status |', '|---|---|---:|---:|---|---|']
        for row in result['readsets']:
            state = row['strand']
            lines.append(f'| {state["assembly_id"]} | {state["readset_id"]} | {row.get("mapped_fraction", 0):.1%} | {state["informative"]} | {state.get("strandedness", "unresolved")} | {state["status"]} |')
            if state.get('reason'):
                lines.extend(['', state['reason'], ''])
        lines.extend(['', '| Assembly | Library | Status | Strand | Detail |', '|---|---|---|---|---|'])
        for row in result['libraries']:
            detail = ', '.join(f'[BAM {i}]({path})' for i, path in enumerate(row['bams'], 1)) if row['status'] == 'ready' else row['reason']
            lines.append(f'| {row["assembly_id"]} | {row["library_id"]} | {row["status"]} | {row.get("strandedness", "unresolved")} | {detail} |')
        lines.extend(['', '[Resolved alignment statistics and BAM inputs](rna_alignment_status.json)', ''])
        Path(args.outdir, 'README.md').write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
