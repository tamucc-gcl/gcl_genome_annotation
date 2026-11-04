#!/usr/bin/env python3
"""RNA alignment contracts and reports. Nextflow stages and transports every file."""
import argparse
import base64
import json
from collections import Counter
from pathlib import Path
from rna_mapping import selected_reads, reference


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def decode(value):
    return json.loads(base64.b64decode(value))


def plan(manifest, prepared, locations):
    if manifest['status'] != 'ready' or prepared['status'] != 'ready':
        raise ValueError('Input preflight and selected RNA preparation must both be ready')
    if locations['strand_policy'].get('unresolved_policy', 'unstranded') not in ('unstranded', 'error'):
        raise ValueError('Invalid unresolved strand policy; use unstranded or error')
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
                         'policy': {**locations['strand_policy'], 'strandedness': current['strandedness']}})
    return {'schema_version': 1, 'assemblies': assemblies, 'references': references, 'jobs': jobs}


def strand(qc, policy):
    requested = policy['strandedness']
    fallback = policy.get('unresolved_policy', 'unstranded')
    if requested not in ('auto', 'forward', 'reverse', 'unstranded'):
        raise ValueError('Invalid strandedness')
    if fallback not in ('unstranded', 'error'):
        raise ValueError('Invalid unresolved strand policy; use unstranded or error')
    minimum, dominant, tolerance = int(policy['min_informative']), float(policy['dominant_fraction']), float(policy['unstranded_tolerance'])
    if minimum < 1 or not .5 < dominant <= 1 or not 0 <= tolerance < dominant - .5:
        raise ValueError('Invalid strand inference thresholds')
    result = {key: qc[key] for key in ('assembly_id', 'readset_id', 'library_id')}
    n = qc['strand_forward'] + qc['strand_reverse']
    fraction = qc['strand_forward'] / n if n else None
    result.update(schema_version=2, status='blocked', alignment_status=qc['status'],
                  requested=requested, requested_strandedness=requested, policy=policy,
                  informative=n, strand_forward=qc['strand_forward'], strand_reverse=qc['strand_reverse'],
                  forward_fraction=fraction, reverse_fraction=qc['strand_reverse'] / n if n else None,
                  method='splice_XS', inferred_strandedness='unresolved',
                  effective_strandedness=None, effective_setting_source=None)
    if qc['status'] != 'ready':
        result['reason'] = qc.get('reason', 'Alignment QC failed')
        result.update(inference_outcome='not_evaluated', inference_reason='Alignment QC failed')
        return result
    if n < minimum:
        outcome, reason = 'insufficient', f'Only {n} informative spliced reads; need {minimum}'
    elif fraction >= dominant:
        outcome, reason = 'confident', 'Forward orientation meets dominant threshold'
        result['inferred_strandedness'] = 'forward'
    elif qc['strand_reverse'] / n >= dominant:
        outcome, reason = 'confident', 'Reverse orientation meets dominant threshold'
        result['inferred_strandedness'] = 'reverse'
    elif abs(fraction - .5) <= tolerance:
        outcome, reason = 'confident', 'Orientation proportions meet unstranded tolerance'
        result['inferred_strandedness'] = 'unstranded'
    else:
        outcome, reason = 'ambiguous', 'Intermediate strand mixture'
    result.update(inference_outcome=outcome, inference_reason=reason)
    if requested != 'auto':
        effective, source = requested, 'explicit'
    elif outcome == 'confident':
        effective, source = result['inferred_strandedness'], 'confident_inference'
    elif fallback == 'unstranded':
        effective, source = 'unstranded', 'conservative_fallback'
        result['warning'] = reason + '; processed without a library-strand assumption'
    else:
        result['reason'] = reason + '; strict strand policy requires confident inference. Supply a known setting or use rna_unresolved_strand_policy=unstranded'
        return result
    # v2 compatibility alias: strandedness means processing, never inferred protocol.
    result.update(status='ready', strandedness=effective, effective_strandedness=effective,
                  effective_setting_source=source)
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
                counts['spliced_xs_reads'] += int(spliced and read.has_tag('XS') and read.get_tag('XS') in ('+', '-'))
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
    result['spliced_xs_reads'] = counts['spliced_xs_reads']
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
    resolved = {(row['assembly_id'], row['readset_id']): dict(row) for row in strands}
    aligned = {(row['assembly_id'], row['readset_id']): row for row in bams}
    libraries = {}
    for job in plan_data['jobs']:
        key = (job['assembly_id'], job['readset_id'])
        if key not in resolved:
            resolved[key] = {**{k: job[k] for k in ('assembly_id', 'readset_id', 'library_id')},
                             'status': 'blocked', 'reason': 'Missing strand analysis',
                             'alignment_status': qc.get(key, {}).get('status', 'missing')}
        libraries.setdefault((job['assembly_id'], job['library_id']), []).append(job['readset_id'])
    rows = []
    for (aid, lid), ids in sorted(libraries.items()):
        states = [resolved.get((aid, rid), {}) for rid in ids]
        paths = [aligned.get((aid, rid), {}) for rid in ids]
        row = dict(assembly_id=aid, library_id=lid, readsets=ids, status='blocked')
        explicit = {state.get('requested', 'auto') for state in states} - {'auto'}
        values = {state.get('effective_strandedness', state.get('strandedness')) for state in states}
        fallback = any(state.get('effective_setting_source') == 'conservative_fallback' for state in states)
        conflict = len(values) != 1 or fallback
        strict = any(state.get('policy', {}).get('unresolved_policy') == 'error' for state in states)
        if len(explicit) > 1:
            row['reason'] = 'Contradictory explicit strand settings within physical library'
        elif not all(state.get('status') == 'ready' for state in states):
            row['reason'] = '; '.join(state.get('reason', 'Missing strand analysis') for state in states if state.get('status') != 'ready')
        elif any(state.get('library_id') != lid or qc.get((aid, rid), {}).get('library_id') != lid
                 for rid, state in zip(ids, states)):
            row['reason'] = 'Readset/library association differs from alignment plan'
        elif any(qc.get((aid, rid), {}).get('status') != 'ready' for rid in ids):
            row['reason'] = 'Alignment QC failed'
        elif not all(path.get('bam') and path.get('bam_index') for path in paths):
            row['reason'] = 'Missing readset BAM/index'
        elif conflict and (explicit or strict):
            row['reason'] = 'Conflicting/unresolved lane strands; review explicit settings or use rna_unresolved_strand_policy=unstranded for automatic libraries'
        else:
            effective = 'unstranded' if conflict else next(iter(values))
            source = ('conservative_fallback' if conflict else 'explicit' if explicit else 'confident_inference')
            row.update(status='ready', strandedness=effective, effective_strandedness=effective,
                       effective_setting_source=source,
                       inferred_strandedness=(next(iter({state.get('inferred_strandedness', 'unresolved') for state in states}))
                           if len({state.get('inferred_strandedness', 'unresolved') for state in states}) == 1 else 'unresolved'),
                       bams=[data['bam'] for data in paths], bam_indices=[data['bam_index'] for data in paths])
            if conflict:
                row['warning'] = 'Automatic lanes disagree or remain unresolved; entire library processed without a library-strand assumption'
            for state in states:
                if conflict:
                    state.update(strandedness=effective, effective_strandedness=effective,
                                 effective_setting_source=source)
                    state['warning'] = state.get('warning', row['warning'])
        rows.append(row)
    assemblies = []
    for row in plan_data['assemblies']:
        relevant = [lib for lib in rows if lib['assembly_id'] == row['assembly_id']]
        assemblies.append({**row, 'status': ('not_requested' if not row['rna_readsets'] else
                           'ready' if relevant and all(lib['status'] == 'ready' for lib in relevant) else 'blocked'),
                           'bams': [path for lib in relevant if lib['status'] == 'ready' for path in lib['bams']]})
    warnings = [dict(assembly_id=state['assembly_id'], library_id=state['library_id'],
                     readset_id=state['readset_id'], warning=state['warning'],
                     informative=state['informative'], strand_forward=state['strand_forward'],
                     strand_reverse=state['strand_reverse'], forward_fraction=state['forward_fraction'],
                     reverse_fraction=state['reverse_fraction'], effective_strandedness=state['effective_strandedness'],
                     effective_setting_source=state['effective_setting_source'],
                     inferred_strandedness=state['inferred_strandedness'], inference_outcome=state['inference_outcome'])
                for state in resolved.values() if state.get('warning')]
    return {'schema_version': 2, 'warnings': warnings, 'has_warnings': bool(warnings),
            'status': 'ready' if all(row['status'] != 'blocked' for row in assemblies) else 'incomplete',
            'assemblies': assemblies, 'libraries': rows,
            'readsets': [{**qc.get((job['assembly_id'], job['readset_id']), {}),
                          'strand': resolved.get((job['assembly_id'], job['readset_id']), {})} for job in plan_data['jobs']]}


def markdown(result):
    lines = ['# RNA alignment', '', f'Status: **{result["status"]}**; warnings: **{len(result["warnings"])}**', '',
             'Ready means operationally usable, not warning-free or biologically validated.',
             'Inference uses unique high-quality spliced primary reads (R1 for pairs); processing settings may use conservative fallbacks.', '',
             '| Assembly | Library | Readset | Alignment | Inferred | Outcome | Effective | Source |',
             '|---|---|---|---|---|---|---|---|']
    for row in result['readsets']:
        state = row['strand']
        lines.append('| ' + ' | '.join(str(state.get(key, 'unresolved')) for key in
            ('assembly_id', 'library_id', 'readset_id', 'alignment_status', 'inferred_strandedness',
             'inference_outcome', 'effective_strandedness', 'effective_setting_source')) + ' |')
        lines.append(f"\n{state['assembly_id']}/{state['library_id']}/{state['readset_id']}: "
                     f"forward={state.get('strand_forward')}, reverse={state.get('strand_reverse')}, informative={state.get('informative')}; "
                     f"fractions={state.get('forward_fraction')}/{state.get('reverse_fraction')}; "
                     + state.get('warning', state.get('reason', state.get('inference_reason', ''))) + '\n')
    for row in result['libraries']:
        lines.append(f"- {row['assembly_id']}/{row['library_id']}: {row['status']}; "
                     + row.get('reason', row.get('warning', row.get('effective_setting_source', ''))))
        for bam in row.get('bams', []):
            lines.append(f'  - [BAM]({bam})')
    lines.extend(['', '[Alignment QC: HISAT2 and samtools](multiqc/multiqc_report.html)',
                  '[Measurements, thresholds, uncertainty and BAM inputs](rna_alignment_status.json)', ''])
    return '\n'.join(lines)


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
        Path(args.outdir, 'README.md').write_text(markdown(result), encoding='utf-8')



if __name__ == '__main__':
    main()
