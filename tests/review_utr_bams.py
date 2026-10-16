#!/usr/bin/env python3
"""Targeted alignment inspection of completed UTR candidates; no filtering."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from evidence_aed import merged
from prepare_protein_evidence import checksum


def library_strand(read, protocol):
    if protocol == 'unstranded':
        return None
    if protocol not in ('forward', 'reverse'):
        raise ValueError('Library strand must be resolved before review')
    reverse = read.is_reverse
    if read.is_paired:
        if read.is_read1 == read.is_read2:
            raise ValueError('Invalid paired-end mate flags')
        if read.is_read2:
            reverse = not reverse
    if protocol == 'reverse':
        reverse = not reverse
    return '-' if reverse else '+'


def cigar_parts(read):
    """Expose matched bases and N gaps with immediate aligned anchor lengths."""
    pos = read.reference_start
    matches, introns = [], []
    cigar = read.cigartuples or []
    for i, (op, length) in enumerate(cigar):
        if op in (0, 7, 8):
            matches.append((pos, pos + length))
        if op == 3:
            left = cigar[i-1][1] if i and cigar[i-1][0] in (0,7,8) else 0
            right = cigar[i+1][1] if i+1 < len(cigar) and cigar[i+1][0] in (0,7,8) else 0
            introns.append((pos, pos + length, min(left, right)))
        if op in (0, 2, 3, 7, 8):
            pos += length
    return matches, introns


def depth_metrics(events, blocks):
    """Read-alignment depth over selected exonic bases, excluding CIGAR N/D."""
    blocks = merged(blocks)
    size = sum(e-s for s,e in blocks)
    points = sorted(set(events) | {x for block in blocks for x in block})
    depth = total = covered = ten = peak = zero_run = longest_zero = 0
    minimum = None
    index = 0
    for i, pos in enumerate(points[:-1]):
        depth += events.get(pos, 0)
        end = points[i+1]
        while index < len(blocks) and blocks[index][1] <= pos:
            index += 1
            zero_run = 0
        if index == len(blocks):
            break
        a,b = blocks[index]
        bases = max(0, min(end,b) - max(pos,a))
        if bases:
            total += bases * depth
            covered += bases if depth else 0
            ten += bases if depth >= 10 else 0
            peak = max(peak,depth)
            minimum = depth if minimum is None else min(minimum,depth)
            zero_run = zero_run + bases if not depth else 0
            longest_zero = max(longest_zero,zero_run)
    return dict(bases=size, mean_depth=total/size if size else None,
        covered_fraction=covered/size if size else None,
        fraction_at_least_10x=ten/size if size else None, maximum_depth=peak,
        minimum_depth=minimum, longest_uncovered_run_bp=longest_zero)


def inspect_bam(bam, case, readset, min_mapq, anchor):
    model = case['accepted']; exons = merged(model['exons'])
    seq, strand = model['seq'], model['strand']
    start, end = exons[0][0], exons[-1][1]
    expected = [(a[1],b[0]) for a,b in zip(exons, exons[1:])]
    original_exons = merged(case['original']['exons'])
    old_junctions = {(a[1],b[0]) for a,b in zip(original_exons,original_exons[1:])}
    events = defaultdict(Counter)
    junctions = defaultdict(lambda: defaultdict(set))
    counts = Counter()
    for read in bam.fetch(seq, start, end):
        counts['fetched_alignments'] += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        counts['primary_mapped_alignments'] += 1
        if read.mapping_quality < min_mapq or not read.has_tag('NH') or read.get_tag('NH') != 1:
            continue
        counts['unique_high_quality_alignments'] += 1
        orientation = library_strand(read,readset['strandedness'])
        group = 'matching_library_strand' if orientation == strand else 'unstranded' if orientation is None else 'opposite_library_strand'
        counts[group] += 1
        if read.has_tag('XS') and orientation and read.get_tag('XS') in ('+','-'):
            counts['library_XS_disagreements'] += read.get_tag('XS') != orientation
        matches, gaps = cigar_parts(read)
        for s,e in matches:
            for key in ('all_strands',group):
                events[key][s] += 1; events[key][e] -= 1
        for s,e,overhang in gaps:
            if overhang >= anchor:
                # One name per gap/library avoids counting both mates as two fragments.
                junctions[(s,e)]['all_strands'].add(read.query_name)
                junctions[(s,e)][group].add(read.query_name)
    segments = []
    for kind, blocks in [('original_exons',case['original']['exons']),
                         ('accepted_extension',case['accepted_extension'])]:
        for i,(s,e) in enumerate(merged(blocks),1):
            segments.append(dict(segment_kind=kind, segment=i,start=s,end=e,
                depth={key:depth_metrics(events[key],[(s,e)]) for key in
                       ('all_strands','matching_library_strand','opposite_library_strand','unstranded')}))
    expected_rows = []
    for s,e in expected:
        support = junctions[(s,e)]
        expected_rows.append(dict(start=s,end=e,new_in_refined_model=(s,e) not in old_junctions,
            supporting_fragments={key:len(support[key]) for key in
                ('all_strands','matching_library_strand','opposite_library_strand','unstranded')}))
    support_key='unstranded' if readset['strandedness']=='unstranded' else 'matching_library_strand'
    new_rows=[j for j in expected_rows if j['new_in_refined_model']]
    return dict(readset_id=readset['readset_id'], protocol=readset['strandedness'],
        counts=dict(counts), segments=segments, expected_junctions=expected_rows,
        new_junction_summary=dict(total=len(new_rows),
            with_matching_support=sum(j['supporting_fragments'][support_key]>0 for j in new_rows),
            without_matching_support=[(j['start'],j['end']) for j in new_rows if not j['supporting_fragments'][support_key]],
            minimum_matching_fragments=min((j['supporting_fragments'][support_key] for j in new_rows),default=None)),
        extension_depth={key:depth_metrics(events[key],case['accepted_extension']) for key in
            ('all_strands','matching_library_strand','opposite_library_strand','unstranded')})


def run(root, aid, review_path, outdir, min_mapq=20, anchor=8):
    import pysam
    root,outdir = Path(root),Path(outdir)
    review = json.loads(Path(review_path).read_text())
    folder = root/'structural_annotation/utr'/aid
    status = json.loads((folder/'utr_status.json').read_text())
    if status['status'] != 'ready' or checksum(folder/f'{aid}.utr.gff3') != status['gff3_sha256']:
        raise ValueError('UTR output differs from accepted status')
    identity = json.loads((root/'filtered_assemblies'/aid/f'{aid}.contigs.json').read_text())
    if identity['reference_identity'] != status['reference_identity']:
        raise ValueError('Reference identity differs from accepted UTR run')
    cases = [c for c in review['selected_cases'] if any(x.startswith('long_') for x in c['selected_for'])]
    # Confirm selected geometry against the accepted file, not an unversioned review JSON.
    from inspect_utr_labels import annotations
    lengths = {c['id']:c['length'] for c in identity['contigs']}
    accepted = annotations(folder/f'{aid}.utr.gff3', lengths)
    original_path=root/'structural_annotation/braker'/aid/f'{aid}.gff3'
    if checksum(original_path) != status['original_model_sha256']:
        raise ValueError('Original annotation differs from accepted UTR run')
    original=annotations(original_path,lengths)
    for c in cases:
        for source, models in [('original',original),('accepted',accepted)]:
            model = models[c['transcript_id']]
            if (model['seq'],model['strand'],sorted(model['exons'])) != (c[source]['seq'],c[source]['strand'],sorted(map(tuple,c[source]['exons']))):
                raise ValueError('Review geometry differs from accepted UTR inputs')
    results = {c['transcript_id']:dict(transcript_id=c['transcript_id'],
        contig=c['accepted']['seq'], strand=c['accepted']['strand'],
        original=c['original'], accepted=c['accepted'], extension=c['accepted_extension'],
        neighboring_genes=c['other_gene_coding_spans_overlapped'], readsets=[]) for c in cases}
    for row in status['readsets']:
        path = root/'rna_alignment'/row['bam']
        index = root/'rna_alignment'/row['bam_index']
        print(f"Inspecting {row['readset_id']} at {len(cases)} selected loci...",flush=True)
        with pysam.AlignmentFile(str(path),'rb',index_filename=str(index)) as bam:
            if list(zip(bam.references,bam.lengths)) != [(c['id'],c['length']) for c in identity['contigs']]:
                raise ValueError('BAM reference dictionary differs from accepted reference')
            if not bam.has_index(): raise ValueError('Indexed BAM required')
            groups = {g['ID']:g for g in bam.header.to_dict().get('RG',[])}
            group = groups.get(row['readset_id'],{})
            if group.get('SM') != row['sample_id'] or group.get('LB') != row['library_id']:
                raise ValueError('BAM read group differs from accepted readset')
            for case in cases:
                results[case['transcript_id']]['readsets'].append(inspect_bam(bam,case,row,min_mapq,anchor))
    outdir.mkdir(parents=True,exist_ok=True)
    result = dict(stage='utr_bam_review',coordinates='0-based half-open',
        model_sha256=status['gff3_sha256'],reference_identity=status['reference_identity'],
        min_mapq=min_mapq,unique_only='NH=1',minimum_immediate_junction_anchor=anchor,
        scope='Local alignment diagnostic, not Portcullis classification or transcript validation; no annotation changes',
        depth_units='aligned read bases; overlapping mates counted separately; no base-quality cutoff; duplicates retained',
        junction_units='distinct query names per junction per readset; independent libraries reported separately',
        cases=list(results.values()))
    (outdir/'bam_review.json').write_text(json.dumps(result,indent=2)+'\n')
    with (outdir/'junction_support.tsv').open('w',newline='') as handle:
        names=['transcript_id','readset_id','contig','strand','start','end','new_in_refined_model',
               'all_strands','matching_library_strand','opposite_library_strand','unstranded']
        w=csv.DictWriter(handle,fieldnames=names,delimiter='\t'); w.writeheader()
        for case in results.values():
            for row in case['readsets']:
                for j in row['expected_junctions']:
                    w.writerow(dict(transcript_id=case['transcript_id'],readset_id=row['readset_id'],
                        contig=case['contig'],strand=case['strand'],start=j['start'],end=j['end'],
                        new_in_refined_model=j['new_in_refined_model'],**j['supporting_fragments']))
    print(f"Completed local BAM review: {outdir}/bam_review.json",flush=True)
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs-root',required=True);p.add_argument('--assembly-id',required=True)
    p.add_argument('--review',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--min-mapq',type=int,default=20);p.add_argument('--min-anchor',type=int,default=8)
    args=p.parse_args()
    if args.min_mapq < 0 or args.min_anchor < 1:p.error('Invalid MAPQ or anchor threshold')
    run(args.inputs_root,args.assembly_id,args.review,args.outdir,args.min_mapq,args.min_anchor)
