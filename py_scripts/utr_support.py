#!/usr/bin/env python3
"""Screen native UTR proposals against the current assembly's RNA alignments."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote, unquote

def table(path, rows, fields):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def gff(path):
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) != 9:
                raise ValueError('Malformed GFF: ' + str(path))
            attrs = {k: unquote(v) for k, v in (x.split('=', 1) for x in c[8].split(';') if '=' in x)}
            yield c, attrs


def merge(blocks):
    result = []
    for s, e in sorted(blocks):
        if s < 0 or e <= s:
            raise ValueError('Invalid interval')
        if result and s <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], e))
        else:
            result.append((s, e))
    return result


def subtract(blocks, excluded):
    parts = merge(blocks)
    for a, b in merge(excluded):
        parts = [(x, y) for s, e in parts for x, y in ((s, min(e, a)), (max(s, b), e)) if y > x]
    return parts


def union_bp(blocks):
    return sum(e - s for s, e in merge(blocks))


def junctions(exons):
    blocks = merge(exons)
    return {(a[1], b[0]) for a, b in zip(blocks, blocks[1:])}


def load_models(path):
    result = {}
    parts = []
    for c, a in gff(path):
        if c[2] in ('mRNA', 'transcript'):
            tid = a['ID']
            if tid in result:
                raise ValueError('Duplicate transcript: ' + tid)
            result[tid] = {'id': tid, 'gene': a['Parent'], 'seq': c[0], 'strand': c[6],
                           'start': int(c[3])-1, 'end': int(c[4]), 'attrs': a,
                           'exons': [], 'cds': [], 'stops': []}
        elif c[2] in ('exon', 'CDS', 'stop_codon'):
            parts.append((c, a))
    for c, a in parts:
        for tid in a.get('Parent', '').split(','):
            if tid not in result:
                raise ValueError('Dangling transcript parent: ' + tid)
            model = result[tid]
            if c[0] != model['seq'] or c[6] != model['strand']:
                raise ValueError('Feature reference/strand mismatch')
            key = {'exon': 'exons', 'CDS': 'cds', 'stop_codon': 'stops'}[c[2]]
            value = (int(c[3])-1, int(c[4]))
            model[key].append(value + (c[7],) if key == 'cds' else value)
    if not result or any(not t['exons'] or not t['cds'] for t in result.values()):
        raise ValueError('Empty/incomplete coding model set')
    return result


def compatible(before, after):
    if set(before) != set(after):
        raise ValueError('Transcript IDs changed')
    for tid, old in before.items():
        new = after[tid]
        if any(old[k] != new[k] for k in ('gene', 'seq', 'strand')) or sorted(old['cds']) != sorted(new['cds']):
            raise ValueError('Coding structure changed: ' + tid)
        if subtract(old['exons'], new['exons']):
            raise ValueError('Existing exon sequence removed: ' + tid)


def overlap_index(before):
    genes = defaultdict(list)
    for t in before.values():
        genes[(t['seq'], t['gene'], t['strand'])].extend((s, e) for s, e, phase in t['cds'])
    bins = defaultdict(list)
    for (seq, gene, strand), blocks in genes.items():
        for s, e in merge(blocks):
            for bucket in range(s//100000, (e-1)//100000+1):
                bins[(seq, bucket)].append((s, e, gene, strand))
    return bins


def geometry(before, after):
    compatible(before, after)
    index = overlap_index(before)
    rows, conflicts = [], []
    for tid, old in before.items():
        new = after[tid]
        added = subtract(subtract(new['exons'], old['exons']), old['stops'])
        lower = min(s for s, e, phase in old['cds'])
        upper = max(e for s, e, phase in old['cds'])
        low = union_bp([(s, min(e, lower)) for s, e in added if s < lower])
        high = union_bp([(max(s, upper), e) for s, e in added if e > upper])
        by_gene = defaultdict(list)
        for s, e in added:
            hits = set(x for bucket in range(s//100000, (e-1)//100000+1)
                       for x in index.get((new['seq'], bucket), []))
            for a, b, gene, strand in hits:
                if gene != old['gene'] and min(e, b) > max(s, a):
                    by_gene[(gene, strand)].append((max(s, a), min(e, b)))
        for (gene, strand), blocks in by_gene.items():
            conflicts.append({'transcript_id': tid, 'other_gene_id': gene,
                              'same_strand': strand == old['strand'], 'overlap_bp': union_bp(blocks)})
        row = {'transcript_id': tid, 'gene_id': old['gene'], 'contig': old['seq'],
               'strand': old['strand'], 'changed': bool(added),
               'old_start_1based': old['start']+1, 'old_end': old['end'],
               'proposed_start_1based': new['start']+1, 'proposed_end': new['end'],
               'added_5prime_spliced_bp': low if old['strand'] == '+' else high,
               'added_3prime_spliced_bp': high if old['strand'] == '+' else low,
               'genomic_left_extension_bp': max(0, old['start']-new['start']),
               'genomic_right_extension_bp': max(0, new['end']-old['end']),
               'added_blocks_0based_halfopen': json.dumps(added),
               'new_junctions': len(junctions(new['exons'])-junctions(old['exons'])),
               'same_strand_cds_overlap_bp': sum(union_bp(v) for (g, st), v in by_gene.items() if st == old['strand']),
               'opposite_strand_cds_overlap_bp': sum(union_bp(v) for (g, st), v in by_gene.items() if st != old['strand']),
               'evidence_id': new['attrs'].get('utr_refine_evidence', '.'),
               'evidence_scope': 'geometry_only', 'decision': 'review' if added else 'unchanged'}
        rows.append(row)
    return rows, conflicts


def effective_read(read, strand, protocol, counters, min_mapq=20):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.mapping_quality < min_mapq:
        counters['low_MAPQ_callbacks']+=1
        return False
    if not read.has_tag('NH'):
        counters['missing_NH'] += 1
        return False
    if read.get_tag('NH') != 1:
        counters['nonunique_callbacks']+=1
        return False
    if protocol != 'unstranded':
        orientation = '-' if read.is_reverse else '+'
        flip = (protocol == 'reverse') == (not read.is_paired or read.is_read1)
        if flip:
            orientation = '+' if orientation == '-' else '-'
        if orientation != strand:
            return False
    counters['accepted_alignment_callbacks'] += 1
    return True


def trial_decision(breadth, gap, minimum, body, missing, ambiguous, conflict,
                   min_breadth=.8, max_gap_bp=100, min_junction_alignments=3):
    if missing:
        return 'review','uniqueness_evidence_incomplete'
    if body>=10 and (breadth==0 or minimum==0) and not ambiguous:
        return 'rejected','trial_no_added_base_or_new_junction_support'
    if breadth>=min_breadth and gap<=max_gap_bp and (minimum is None or minimum>=min_junction_alignments) and not conflict:
        return 'accepted','trial_support_thresholds_met'
    return 'review','trial_weak_ambiguous_support_or_same_strand_conflict'


def support(identity, metadata, bams, indices, before, after, rows, rules):
    import numpy as np
    import pysam

    anchor = rules['anchor_bp']
    dictionary = [(row['id'], row['length']) for row in identity['contigs']]
    evidence = list(zip(metadata, bams, indices, strict=True))
    for _, bam, index in evidence:
        with pysam.AlignmentFile(str(bam), 'rb', index_filename=str(index)) as handle:
            if list(zip(handle.references, handle.lengths)) != dictionary:
                raise ValueError('BAM reference dictionary differs from accepted reference')
    details = []
    for row in rows:
        tid = row['transcript_id']
        if not row['changed']:
            continue
        old, new = before[tid], after[tid]
        blocks = subtract(subtract(new['exons'], old['exons']), old['stops'])
        totals = [np.zeros(e-s, dtype=np.uint64) for s, e in blocks]
        junction_totals = Counter()
        body_total = 0
        missing = 0
        ambiguous = 0
        new_junctions = junctions(new['exons'])-junctions(old['exons'])
        boundaries=set()
        old_blocks=old['exons']+old['stops']
        for point in (min(s for s,e in old_blocks),max(e for s,e in old_blocks)):
            if any(s<point<e for s,e in merge(new['exons'])) and any(s<=point<=e for s,e in blocks):
                boundaries.add(point)
        library_depth={};library_junctions=defaultdict(Counter);library_boundaries=defaultdict(Counter)
        for meta, bam, index in evidence:
            library=meta['library_id']
            if library not in library_depth:library_depth[library]=[np.zeros(e-s,dtype=np.uint64) for s,e in blocks]
            protocol = meta.get('effective_strandedness', meta['strandedness'])
            if protocol not in ('forward', 'reverse', 'unstranded'):
                raise ValueError('Unresolved effective protocol')
            stats = Counter()
            callback = lambda read: effective_read(read, old['strand'], protocol, stats, rules['min_mapq'])
            counts = Counter()
            with pysam.AlignmentFile(str(bam), 'rb', index_filename=str(index)) as handle:
                for i, (s, e) in enumerate(blocks):
                    arrays = handle.count_coverage(old['seq'], s, e, quality_threshold=0, read_callback=callback)
                    depth = np.sum(np.asarray(arrays, dtype=np.uint64), axis=0)
                    totals[i] += depth
                    library_depth[library][i]+=depth
                    counts['covered_bases'] += int(np.count_nonzero(depth))
                s, e, phase = sorted(old['cds'])[0]
                counts['coding_alignments'] = sum(callback(r) for r in handle.fetch(old['seq'], s, min(e, s+500)))
                for donor, acceptor in new_junctions:
                    for alignment in handle.fetch(old['seq'], max(0, donor-1), donor+1):
                        if not callback(alignment):
                            continue
                        pos = alignment.reference_start
                        for j, (op, length) in enumerate(alignment.cigartuples or []):
                            if op == 3 and (pos, pos+length) == (donor, acceptor):
                                cigar = alignment.cigartuples
                                if j > 0 and j+1 < len(cigar) and cigar[j-1][0] in (0, 7, 8) and cigar[j+1][0] in (0, 7, 8) and min(cigar[j-1][1], cigar[j+1][1]) >= anchor:
                                    junction_totals[(donor, acceptor)] += 1
                                    library_junctions[library][(donor,acceptor)]+=1
                                    counts['new_junction_alignment_occurrences'] += 1
                            if op in (0, 2, 3, 7, 8):
                                pos += length
                for point in boundaries:
                    for alignment in handle.fetch(old['seq'],max(0,point-anchor),point+anchor):
                        if callback(alignment) and any(s<=point-anchor and e>=point+anchor for s,e in alignment.get_blocks()):
                            library_boundaries[library][point]+=1
            missing += stats['missing_NH']
            ambiguous += stats['low_MAPQ_callbacks']+stats['nonunique_callbacks']
            body_total += counts['coding_alignments']
            details.append({'transcript_id': tid, 'readset_id': meta['readset_id'],
                            'library_id': meta['library_id'], 'effective_strandedness': protocol,
                            'missing_NH_callbacks': stats['missing_NH'],
                            'ambiguous_mapping_callbacks':stats['low_MAPQ_callbacks']+stats['nonunique_callbacks'],**counts})
        denominator = sum(len(x) for x in totals)
        breadth = sum(int(np.count_nonzero(x)) for x in totals)/denominator if denominator else 0
        longest_gap = 0
        for depth in totals:
            gap = 0
            for value in depth:
                gap = gap+1 if value == 0 else 0
                longest_gap = max(longest_gap, gap)
        minimum = min((junction_totals[j] for j in new_junctions), default=None)
        row.update(evidence_scope='current_run_RNA_alignments', covered_fraction=breadth,
                   longest_uncovered_exonic_gap_bp=longest_gap, coding_alignment_occurrences=body_total,
                   min_new_junction_alignments=minimum if minimum is not None else '.', missing_NH_callbacks=missing,
                   ambiguous_mapping_callbacks=ambiguous)
        library_rows=[]
        for library,depths in library_depth.items():
            library_rows.append({'library_id':library,'covered_fraction':sum(int(np.count_nonzero(x)) for x in depths)/denominator,
                'all_new_junctions_supported':all(library_junctions[library][j]>=rules['min_junction_alignments'] for j in new_junctions),
                'all_boundary_links_supported':all(library_boundaries[library][p]>=rules['min_boundary_alignments'] for p in boundaries),
                'min_new_junction_alignments':min((library_junctions[library][j] for j in new_junctions),default=None),
                'boundary_alignment_counts':{str(p):library_boundaries[library][p] for p in sorted(boundaries)}})
        row['library_support_json']=json.dumps(library_rows,sort_keys=True)
        row['coherent_library_count']=sum(x['covered_fraction']>=rules['min_breadth'] and x['all_new_junctions_supported'] and x['all_boundary_links_supported'] for x in library_rows)
        row['best_library_covered_fraction']=max(x['covered_fraction'] for x in library_rows)
        row['boundary_links_assessed']=len(boundaries)
        row['terminal_boundary_validation']='not_established_by_short_read_support'
        # Trial policy only: missing uniqueness tags cannot become evidence of absence.
        decision,reason=trial_decision(breadth,longest_gap,minimum,body_total,missing,ambiguous,row['same_strand_cds_overlap_bp'],
            rules['min_breadth'],rules['max_gap_bp'],rules['min_junction_alignments'])
        row.update(decision=decision,decision_reason=reason)
        flags=[]
        if missing:flags.append('uniqueness_evidence_incomplete')
        if breadth<rules['min_breadth']:flags.append('covered_fraction_below_trial_minimum')
        if longest_gap>rules['max_gap_bp']:flags.append('uncovered_gap_above_trial_maximum')
        if minimum is not None and minimum<rules['min_junction_alignments']:flags.append('new_junction_support_below_trial_minimum')
        if row['same_strand_cds_overlap_bp']:flags.append('same_strand_CDS_assignment_conflict')
        if ambiguous and (breadth==0 or minimum==0):flags.append('mapping_ambiguity_prevents_absence_rejection')
        if not row['coherent_library_count']:flags.append('no_coherent_library')
        row['trial_rule_flags']=';'.join(flags) or '.'
        row['coherence_trial_decision']='review' if decision=='accepted' and row['coherent_library_count']==0 else decision
    return details


POLICY = dict(min_breadth=.8, max_gap_bp=100, min_junction_alignments=3,
              min_boundary_alignments=3, anchor_bp=8, min_mapq=20)


def validate_policy(values):
    if set(values)-set(POLICY):raise ValueError('Unknown UTR support setting')
    result=dict(POLICY,**values)
    for key,value in result.items():
        if not isinstance(value,(float,int)) or isinstance(value,bool) or not math.isfinite(value):raise ValueError('Non-numeric UTR setting')
        if key=='min_breadth':
            if not 0<value<=1:raise ValueError('UTR breadth must be >0 and <=1')
        elif value!=int(value) or value<(0 if key in ('max_gap_bp','min_mapq') else 1):
            raise ValueError('Invalid UTR count/length setting: '+key)
        else:result[key]=int(value)
    return result


def accepted(row):
    """Support proxy plus coding assignment; no full-length validation claim."""
    return (row['changed'] and row.get('decision') == 'accepted'
            and int(row.get('coherent_library_count', 0)) > 0
            and int(row['same_strand_cds_overlap_bp']) == 0)


def project(original, proposal, chosen, output):
    """Project whole transcript subtrees; preserve codons and reset gene bounds."""
    models = load_models(original)
    compatible(models, load_models(proposal))
    if set(chosen) != set(models) or set(chosen.values()) - {'original', 'new'}:
        raise ValueError('Every transcript must select the original or current proposal')
    genes = {}; tx = {}; children = defaultdict(list)
    for label, path in (('original', original), ('new', proposal)):
        for c, a in gff(path):
            if c[2] == 'gene':
                if label == 'original': genes[a['ID']] = (list(c), a)
            elif c[2] in ('mRNA', 'transcript'):
                if chosen[a['ID']] == label: tx[a['ID']] = (list(c), dict(a))
            elif 'Parent' in a:
                for tid in a['Parent'].split(','):
                    if tid not in models or chosen[tid] != label: continue
                    attrs = dict(a, Parent=tid, ID='trial.'+tid+'.feature'+str(len(children[tid])+1))
                    children[tid].append((list(c), attrs))
    if set(tx) != set(models): raise ValueError('Candidate projection lost transcripts')
    by_gene = defaultdict(list)
    for tid, m in models.items(): by_gene[m['gene']].append(tid)
    def line(c, a):
        return '\t'.join(c[:8]+[';'.join(k+'='+quote(str(v),safe=',:._-') for k,v in a.items())])+'\n'
    with Path(output).open('x') as h:
        h.write('##gff-version 3\n')
        for gid, tids in by_gene.items():
            c,a=genes[gid]; c[3]=str(min(int(tx[t][0][3]) for t in tids)); c[4]=str(max(int(tx[t][0][4]) for t in tids))
            h.write(line(c,a))
            for tid in tids:
                c,a=tx[tid]; h.write(line(c,a))
                for c,a in children[tid]: h.write(line(c,a))
    after=load_models(output);compatible(models,after)
    parts=defaultdict(lambda:defaultdict(list))
    for c,a in gff(output):
        if c[2] in ('five_prime_UTR','three_prime_UTR'):
            for tid in a['Parent'].split(','):parts[tid][c[2]].append((int(c[3])-1,int(c[4])))
    for tid in models:
        # Compare codons independently of CDS signatures.
        if sorted(models[tid]['stops']) != sorted(after[tid]['stops']):
            raise ValueError('Projection changed stop codon: '+tid)
        if chosen[tid]!='original':
            m=after[tid];lower=min(a for a,b,phase in m['cds']);upper=max(b for a,b,phase in m['cds'])
            low=[(a,min(b,lower)) for a,b in m['exons'] if a<lower]
            high=[(max(a,upper),b) for a,b in m['exons'] if b>upper]
            expected={'five_prime_UTR':low if m['strand']=='+' else high,
                      'three_prime_UTR':high if m['strand']=='+' else low}
            for kind,blocks in expected.items():
                if merge(parts[tid][kind])!=subtract(blocks,m['stops']):
                    raise ValueError('Candidate UTR labels/orientation disagree with exons: '+tid)
    def codons(path):
        result=defaultdict(list)
        for c,a in gff(path):
            if c[2] in ('start_codon','stop_codon'):
                for tid in a['Parent'].split(','): result[tid].append(tuple(c[i] for i in (0,2,3,4,6,7)))
        return {t:sorted(v) for t,v in result.items()}
    if codons(original)!=codons(output):raise ValueError('Projection changed original codons')



def screen(original, proposed, identity, inputs, bams, indices, native_audit, output, policy):
    inputs = json.loads(Path(inputs).read_text())
    identity = json.loads(Path(identity).read_text())
    metadata = inputs['readsets']
    if not metadata:
        raise ValueError('UTR screening requires the current run RNA readsets')
    before, after = load_models(original), load_models(proposed)
    rows, conflicts = geometry(before, after)
    detail = support(identity, metadata, bams, indices, before, after, rows, validate_policy(policy))
    chosen = {row['transcript_id']: 'new' if accepted(row) else 'original' for row in rows}
    output = Path(output)
    output.mkdir()
    project(original, proposed, chosen, output/'screened.gff3')
    audit = json.loads(Path(native_audit).read_text())
    audit['retained_original_transcripts'] = [tid for tid, label in chosen.items() if label == 'original']
    audit['support_withheld'] = {
        row['transcript_id']: row['trial_rule_flags'] for row in rows
        if row['changed'] and chosen[row['transcript_id']] == 'original'
    }
    audit['support_policy'] = validate_policy(policy)
    (output/'screened.rejections.json').write_text(json.dumps(audit,indent=2)+'\n')
    for row in rows:
        row['selected_model'] = chosen[row['transcript_id']]
    table(output/'utr_support.tsv', rows, list(dict.fromkeys(k for row in rows for k in row)))
    table(output/'utr_proposal_conflicts.tsv', conflicts,
          ['transcript_id','other_gene_id','same_strand','overlap_bp'])
    table(output/'utr_readsets.tsv', detail, list(dict.fromkeys(k for row in detail for k in row)) or ['transcript_id'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('original','proposed','identity','inputs','native-audit','outdir','policy'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--bams', nargs='+', required=True)
    parser.add_argument('--indices', nargs='+', required=True)
    args = parser.parse_args()
    screen(args.original,args.proposed,args.identity,args.inputs,args.bams,args.indices,
           args.native_audit,args.outdir,json.loads(args.policy))
