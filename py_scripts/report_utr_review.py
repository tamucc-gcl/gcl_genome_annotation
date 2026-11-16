#!/usr/bin/env python3
"""Report strand evidence for released UTR extensions."""

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import re
from urllib.parse import unquote

from reference_repeats import sha
from utr_support import effective_read, junctions, load_models, merge, overlap_index, table


def read_table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_models(path):
    models = defaultdict(list)
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t', 8)
            if len(fields) != 9 or fields[2] != 'exon' or fields[6] not in ('+', '-'):
                continue
            gtf = re.search(r'(?:^|;)\s*transcript_id\s+"([^"]+)"', fields[8])
            parents = [gtf[1]] if gtf else next(
                ([unquote(parent) for parent in item[7:].split(',')]
                 for item in fields[8].split(';') if item.startswith('Parent=')), [])
            for transcript in parents:
                models[fields[0], fields[6], transcript].append((int(fields[3])-1, int(fields[4])))
    return {key: merge(blocks) for key, blocks in models.items()}


def mirror_matches(rows, models):
    by_chain = defaultdict(list)
    for (contig, strand, name), blocks in models.items():
        chain = tuple(sorted(junctions(blocks)))
        for bucket in range(blocks[0][0]//100000, (blocks[-1][1]-1)//100000+1):
            by_chain[contig, strand, chain, bucket].append(name)
    matches, missing = defaultdict(list), set()
    for tid, row in rows.items():
        contig, strand, evidence = (row[key] for key in ('contig', 'strand', 'evidence_id'))
        blocks = models.get((contig, strand, evidence))
        if blocks is None:
            missing.add(tid)
            continue
        chain = tuple(sorted(junctions(blocks)))
        opposite = '-' if strand == '+' else '+'
        candidates = {name for bucket in range(blocks[0][0]//100000, (blocks[-1][1]-1)//100000+1)
                      for name in by_chain.get((contig, opposite, chain, bucket), [])}
        for name in sorted(candidates):
            other = models[contig, opposite, name]
            if blocks[0][0] >= other[-1][1] or other[0][0] >= blocks[-1][1]:
                continue
            # Matching intron coordinates leave only the outer exon ends to compare.
            lengths = sum(b-a for a, b in blocks), sum(b-a for a, b in other)
            shared = sum(min(b,d)-max(a,c) for (a,b), (c,d) in zip(blocks, other))
            if shared * 100 >= 95 * max(lengths):
                matches[tid].append((name, 'exact' if blocks == other else 'near'))
    return matches, missing


def cds_overlaps(index, row):
    hits = defaultdict(list)
    for start, end in json.loads(row['added_blocks_0based_halfopen']):
        candidates = {hit for bucket in range(start//100000, (end-1)//100000+1)
                      for hit in index.get((row['contig'], bucket), [])}
        for a, b, gene, strand in candidates:
            if gene != row['gene_id'] and strand != row['strand'] and max(a,start) < min(b,end):
                hits[gene].append((max(a,start), min(b,end)))
    intervals = [(gene, a, b) for gene in sorted(hits) for a, b in merge(hits[gene])]
    if sum(b-a for _, a, b in intervals) != int(row['opposite_strand_cds_overlap_bp']):
        raise ValueError('Reported opposite-strand CDS overlap differs for '+row['transcript_id'])
    return intervals


def junction_counts(bam, contig, junction, strand, protocol, policy):
    start, end = junction
    counts, unused = Counter(), Counter()
    anchor = policy['anchor_bp']
    for read in bam.fetch(contig, max(0, start-1), start+1):
        if not effective_read(read, strand, 'unstranded', unused, policy['min_mapq']):
            continue
        pos = read.reference_start
        cigar = read.cigartuples or []
        for i, (op, length) in enumerate(cigar):
            if op == 3 and (pos, pos+length) == (start, end) and 0 < i < len(cigar)-1:
                before, after = cigar[i-1], cigar[i+1]
                if before[0] in (0, 7, 8) and after[0] in (0, 7, 8) and min(before[1], after[1]) >= anchor:
                    if protocol == 'unstranded':
                        counts['unstranded'] += 1
                    elif effective_read(read, strand, protocol, unused, policy['min_mapq']):
                        counts['proposed'] += 1
                    else:
                        counts['opposite'] += 1
            if op in (0, 2, 3, 7, 8):
                pos += length
    return counts


def overlap_depth(bam, contig, start, end, strand, protocol, mapq):
    import numpy as np

    counters = Counter()
    arrays = bam.count_coverage(contig, start, end, quality_threshold=0,
        read_callback=lambda read: effective_read(read, strand, protocol, counters, mapq))
    depth = np.sum(arrays, axis=0, dtype=np.uint64)
    return int(np.count_nonzero(depth)), int(depth.sum())


def protocol_of(meta):
    return meta.get('effective_strandedness') or meta['strandedness']


def junction_flags(counts, stranded, minimum):
    if not stranded:
        return {'stranded_rna_unavailable_for_new_junction'}
    proposed = [c['proposed'] for c in counts.values()]
    total = sum(proposed)
    if not total:
        return {'new_junction_supported_only_by_unstranded_reads'}
    if max(proposed) < minimum:
        if total >= minimum:
            return {'weak_stranded_junction_support_across_libraries'}
        return {'weak_stranded_junction_support'}
    return set()


def review(original, screened, support, clustered, release, status_path, bams, indices, out):
    status = json.loads(status_path.read_text(encoding='utf-8'))
    release_status = json.loads((release/'release_status.json').read_text(encoding='utf-8'))
    if (release_status['status'] != 'ready' or release_status['reference_identity'] != status['reference_identity']
            or sha(original) != status['original_model_sha256'] or sha(screened) != status['gff3_sha256']):
        raise ValueError('UTR review inputs differ from the accepted annotation')
    id_path = 'id_maps/transcript_id_map.tsv'
    gene_path = 'id_maps/gene_id_map.tsv'
    reference_path = 'provenance/reference.json'
    for path in (id_path, gene_path, reference_path):
        if sha(release/path) != release_status['outputs'][path]['sha256']:
            raise ValueError('Release ID map changed: '+path)
    ids = {r['transcript_id']: r for r in read_table(release/id_path)}
    genes = {r['gene_id']: r['final_gene_id'] for r in read_table(release/gene_path)}
    rows = {r['transcript_id']: r for r in read_table(support)
            if r['selected_model'] == 'new' and ids[r['transcript_id']]['decision'] == 'retained'}
    references = status['readsets']
    dictionary = [(row['id'], row['length']) for row in
                  json.loads((release/reference_path).read_text(encoding='utf-8'))['contigs']]
    protocols = {r['library_id']: set() for r in references}
    for r in references:
        protocols[r['library_id']].add(protocol_of(r))
    stranded = {lib for lib, values in protocols.items() if values & {'forward', 'reverse'}}
    mixed = {lib for lib in stranded if 'unstranded' in protocols[lib]}
    pure_stranded = stranded - mixed
    policy = status['utr_support_policy']
    minimum = policy['min_junction_alignments']
    mirrors, missing = mirror_matches(rows, read_models(clustered))
    before, after = load_models(original), load_models(screened)
    cds = overlap_index(before)
    flags, detail, coverage = [], [], []
    import pysam
    with ExitStack() as stack:
        alignments = [(meta, stack.enter_context(pysam.AlignmentFile(str(bam), 'rb', index_filename=str(index))))
                      for meta, bam, index in zip(references, bams, indices, strict=True)]
        for meta, bam in alignments:
            if list(zip(bam.references, bam.lengths)) != dictionary:
                raise ValueError('UTR review BAM reference dictionary differs for '+meta['readset_id'])
        for tid, row in rows.items():
            libraries = json.loads(row['library_support_json'])
            passing = [r['library_id'] for r in libraries if r['covered_fraction'] >= policy['min_breadth']
                       and r['all_new_junctions_supported'] and r['all_boundary_links_supported']]
            reasons = set()
            if tid in mirrors:
                reasons.add('opposite_strand_rna_model')
                if not set(passing) & pure_stranded:
                    reasons.add('opposite_strand_rna_model_without_stranded_library_pass')
            if int(row['opposite_strand_cds_overlap_bp']) and not set(passing) & pure_stranded:
                reasons.add('opposite_cds_without_stranded_library_pass')
            junction_rows = []
            contig, strand = row['contig'], row['strand']
            opposite = '-' if strand == '+' else '+'
            junction_review = not any(r['library_id'] in pure_stranded and r['all_new_junctions_supported']
                                      for r in libraries)
            if int(row['new_junctions']) and (reasons or junction_review):
                new = sorted(junctions(after[tid]['exons']) - junctions(before[tid]['exons']))
                if len(new) != int(row['new_junctions']):
                    raise ValueError('New junction count differs for '+tid)
                joint_support = stranded.copy()
                problems = set()
                observed_minimum = {}
                for start, end in new:
                    counts = defaultdict(Counter)
                    for meta, bam in alignments:
                        protocol = protocol_of(meta)
                        counts[meta['library_id']].update(junction_counts(bam,contig,(start,end),strand,protocol,policy))
                    problems.update(junction_flags(counts, stranded, minimum))
                    joint_support.intersection_update(lib for lib, c in counts.items() if c['proposed'] >= minimum)
                    for library, c in sorted(counts.items()):
                        total = c['proposed'] + c['unstranded']
                        observed_minimum[library] = min(observed_minimum.get(library, total), total)
                        junction_rows.append([start,end,library,','.join(sorted(protocols[library])),
                            c['proposed'] if library in stranded else '.',
                            c['opposite'] if library in stranded else '.',
                            c['unstranded'] if 'unstranded' in protocols[library] else '.'])
                expected_minimum = {r['library_id']: r['min_new_junction_alignments'] for r in libraries}
                if observed_minimum != expected_minimum:
                    raise ValueError('Review junction counts disagree with the UTR screen for '+tid)
                if not problems and not joint_support:
                    problems.add('new_junction_support_split_across_stranded_libraries')
                reasons.update(problems)
            if not reasons:
                continue
            release_id = ids[tid]['final_transcript_id']
            flags.append(dict(release_transcript_id=release_id, original_transcript_id=tid,
                contig=contig,strand=strand,flags=','.join(sorted(reasons)),evidence_id=row['evidence_id'],
                added_5prime_spliced_bp=row['added_5prime_spliced_bp'],added_3prime_spliced_bp=row['added_3prime_spliced_bp'],
                new_junctions=row['new_junctions'],opposite_strand_cds_overlap_bp=row['opposite_strand_cds_overlap_bp'],
                stranded_libraries_passing=','.join(lib for lib in passing if lib in pure_stranded),
                unstranded_libraries_passing=','.join(lib for lib in passing if lib not in stranded),
                mixed_libraries_passing=','.join(lib for lib in passing if lib in mixed),
                opposite_strand_rna_models=','.join(f'{name}:{kind}' for name, kind in sorted(mirrors.get(tid, [])))))
            detail.extend([release_id,tid,contig,strand,*r] for r in junction_rows)
            intervals = cds_overlaps(cds, row)
            for gene, start, end in intervals:
                for meta, bam in alignments:
                    protocol = protocol_of(meta)
                    for label, direction in ([('unstranded', strand)] if protocol == 'unstranded'
                                             else [('proposed', strand), ('opposite', opposite)]):
                        covered, depth = overlap_depth(bam,contig,start,end,direction,protocol,policy['min_mapq'])
                        coverage.append([release_id,tid,gene,genes.get(gene,'.'),contig,start,end,
                            meta['library_id'],meta['readset_id'],protocol,label,
                            direction if label != 'unstranded' else '.',covered,depth])

    flag_fields = ('release_transcript_id original_transcript_id contig strand flags evidence_id '
                   'added_5prime_spliced_bp added_3prime_spliced_bp new_junctions opposite_strand_cds_overlap_bp '
                   'stranded_libraries_passing unstranded_libraries_passing mixed_libraries_passing '
                   'opposite_strand_rna_models').split()
    junction_fields = ('release_transcript_id original_transcript_id contig proposed_strand intron_start_0based '
                       'intron_end_exclusive library protocols proposed_strand_alignments opposite_strand_alignments '
                       'unstranded_alignments').split()
    coverage_fields = ('release_transcript_id original_transcript_id opposite_original_gene_id opposite_release_gene_id '
                       'contig start_0based end_exclusive library readset protocol coverage_strand rna_strand covered_bp '
                       'aligned_base_occurrences').split()
    out.mkdir()
    table(out/'utr_review.tsv', sorted(flags, key=lambda r: r['release_transcript_id']), flag_fields)
    table(out/'utr_junction_support.tsv', (dict(zip(junction_fields,r)) for r in detail), junction_fields)
    table(out/'utr_cds_overlap_coverage.tsv', (dict(zip(coverage_fields,r)) for r in coverage), coverage_fields)
    exact = sum(any(kind == 'exact' for _,kind in hits) for hits in mirrors.values())
    summary = dict(released_extensions=len(rows), flagged_transcripts=len(flags), mirror_exact=exact,
                   mirror_near_only=len(mirrors)-exact, mirror_evidence_models_missing=len(missing),
                   stranded_library_count=len(stranded), flags=dict(Counter(reason for row in flags
                   for reason in row['flags'].split(','))), reference_identity=status['reference_identity'],
                   release_status_sha256=sha(release/'release_status.json'),
                   source_sha256=dict(original_model=status['original_model_sha256'], screened_model=status['gff3_sha256'],
                                      utr_status=sha(status_path), utr_support=sha(support), clustered_rna=sha(clustered)))
    (out/'utr_review_status.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('original', 'screened', 'support', 'clustered', 'release', 'status', 'outdir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--bams', nargs='+', type=Path, required=True)
    parser.add_argument('--indices', nargs='+', type=Path, required=True)
    args = parser.parse_args()
    review(args.original,args.screened,args.support,args.clustered,args.release,args.status,args.bams,args.indices,args.outdir)


if __name__ == '__main__':
    main()
