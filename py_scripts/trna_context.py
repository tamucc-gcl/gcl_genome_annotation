#!/usr/bin/env python3
"""Report RNA and repeat context for selected tRNA loci."""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

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


def repeat_class(label):
    import re
    if re.search(r'SINE[/_-]tRNA|tRNA[-_ ]derived',label,re.I):return 'trna_derived_TE'
    if re.search(r'\bSINE\b|\bLINE\b|\bLTR\b|DNA/|RC/|Helitron',label,re.I):return 'other_TE'
    if re.search(r'\btRNA\b|\brRNA\b|\bsnRNA\b|\bscRNA\b|\bRNA\b',label,re.I):return 'RNA_label'
    if re.search(r'Satellite|satDNA',label,re.I):return 'satellite'
    return 'other_or_unclassified'

def report(selected, native, decisions, repeats, identity, out):
    with Path(decisions).open() as handle:
        categories = {row['feature_id']:row['native_category'] for row in csv.DictReader(handle,delimiter='\t')}
    native_rows = {}
    with Path(native).open() as handle:
        for line in handle:
            cols = [v.strip() for v in line.rstrip().split('\t')]
            if len(cols) < 13 or not all(cols[i].isdigit() for i in (1,2,3)):
                continue
            start,end = int(cols[2]),int(cols[3])
            key = (cols[0],min(start,end)-1,max(start,end),'+' if start<end else '-')
            if key in native_rows:
                raise ValueError('Duplicate native tRNA locus')
            native_rows[key] = cols

    lengths = {r['id']:r['length'] for r in json.loads(Path(identity).read_text())['contigs']}
    calls = []
    for cols,attrs in gff(selected):
        if cols[2] != 'tRNA':
            continue
        seq,start,end,strand = cols[0],int(cols[3])-1,int(cols[4]),cols[6]
        row = native_rows.get((seq,start,end,strand))
        if row is None:
            raise ValueError('Selected tRNA lacks native score row: '+attrs['ID'])
        calls.append(dict(feature_id=attrs['ID'],contig=seq,start=start,end=end,strand=strand,
            isotype=row[4],anticodon=row[5],domain_score=row[8],structure_score=row[10],
            isotype_model=row[11],isotype_score=row[12],native_classification=categories.get(attrs['ID'],'.'),
            contig_edge_distance=min(start,lengths[seq]-end)))
    calls.sort(key=lambda r:(r['contig'],r['start'],r['end'],r['feature_id']))
    groups = []
    bins = defaultdict(set)
    copies = Counter((row['isotype'],row['anticodon']) for row in calls)
    for i,row in enumerate(calls):
        if not groups or groups[-1][0]!=row['contig'] or row['start']-groups[-1][1]>1000:
            groups.append([row['contig'],row['end'],[]])
        group = groups[-1]
        group[1] = max(group[1],row['end'])
        group[2].append(i)
        for b in range(row['start']//100000,(row['end']-1)//100000+1):
            bins[(row['contig'],b)].add(i)

    blocks = defaultdict(lambda: defaultdict(list))
    labels = defaultdict(set)
    if bins and repeats:
        for cols,attrs in gff(repeats):
            start,end = int(cols[3])-1,int(cols[4])
            ids = set().union(*(bins.get((cols[0],b),set()) for b in range(start//100000,(end-1)//100000+1)))
            if not ids:
                continue
            label = ' '.join(attrs.get(k,'') for k in ('Name','Target','class','classification','family','repeat_class'))
            kind = repeat_class(label)
            for i in ids:
                row = calls[i]
                lo,hi = max(start,row['start']),min(end,row['end'])
                if lo < hi:
                    blocks[i]['repeat'].append((lo,hi))
                    blocks[i][kind].append((lo,hi))
                    labels[i].add(label.strip())

    for j,(_,_,members) in enumerate(groups,1):
        for i in members:
            calls[i].update(cluster_id=f'cluster_{j:05d}',cluster_loci=len(members))
    classes = ('repeat','RNA_label','trna_derived_TE','other_TE','satellite','other_or_unclassified')
    for i,row in enumerate(calls):
        row.update(copy_count=copies[(row['isotype'],row['anticodon'])],repeat_labels=';'.join(sorted(labels[i])) or '.')
        for kind in classes:
            bp = union_bp(blocks[i][kind])
            row[kind+'_overlap_bp'] = bp if repeats else '.'
            row[kind+'_overlap_fraction'] = bp/(row['end']-row['start']) if repeats else '.'
        flags = []
        if blocks[i]['repeat']: flags.append('repeat_overlap_context')
        if blocks[i]['trna_derived_TE']: flags.append('explicit_trna_derived_repeat_overlap')
        if row['cluster_loci']>1: flags.append('local_cluster_context')
        if row['copy_count']>40: flags.append('copy_context_gt40')
        if row['contig_edge_distance']<100: flags.append('near_contig_edge')
        row['context_flags'] = ';'.join(flags) or '.'
    fields = ['feature_id','contig','start','end','strand','isotype','anticodon','domain_score','structure_score',
        'isotype_model','isotype_score','native_classification','contig_edge_distance','cluster_id','cluster_loci',
        'copy_count','repeat_labels']+[kind+suffix for kind in classes for suffix in ('_overlap_bp','_overlap_fraction')]+['context_flags']
    with Path(out).open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fields,delimiter='\t')
        writer.writeheader()
        writer.writerows(calls)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('selected','native','decisions','identity','out'):
        p.add_argument('--'+name,required=True)
    p.add_argument('--repeats')
    a = p.parse_args()
    report(a.selected,a.native,a.decisions,a.repeats,a.identity,a.out)
