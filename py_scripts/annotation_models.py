"""Shared GFF/GTF model, RNA-chain and interval validation (half-open coordinates)."""
from pathlib import Path
import re
from urllib.parse import unquote

def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError('Invalid half-open interval')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result

def overlap_length(left, right):
    left, right = merged(left), merged(right)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total

def features(path, lengths, gtf=False):
    with Path(path).open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            if line.startswith('>'):
                raise ValueError('Embedded FASTA is not supported in evidence GFF')
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected 9 fields')
            seq, _, kind, begin, finish, _, strand, phase, text = fields
            start, end = int(begin) - 1, int(finish)
            if seq not in lengths or not 0 <= start < end <= lengths[seq]:
                raise ValueError(f'{path}:{number}: feature outside accepted reference')
            if strand not in (('+', '-', '.') if gtf else ('+', '-')):
                raise ValueError(f'{path}:{number}: unresolved feature strand')
            attrs = (dict(re.findall(r'(\w+)\s+"([^"\n]*)"', text)) if gtf else
                     {k: unquote(v) for k, v in (x.split('=', 1) for x in text.split(';') if '=' in x)})
            yield dict(seq=seq, kind=kind, start=start, end=end, strand=strand, phase=phase, attrs=attrs)

def models(path, lengths):
    genes, transcripts, parts = {}, {}, []
    for row in features(path, lengths):
        kind, attrs = row['kind'], row['attrs']
        identifier = attrs.get('ID')
        parents = attrs.get('Parent', '').split(',')
        if kind == 'gene':
            if not identifier or identifier in genes:
                raise ValueError('Duplicate/missing model gene ID')
            genes[identifier] = row
        elif kind in ('mRNA', 'transcript'):
            if not identifier or identifier in transcripts or len(parents) != 1 or not parents[0]:
                raise ValueError('Model transcript requires unique ID and explicit gene parent')
            transcripts[identifier] = row | dict(gene_id=parents[0], exons=[], cds=[])
        elif kind in ('exon', 'CDS'):
            if not parents[0]:
                raise ValueError('Model block missing Parent')
            parts.append((row, parents))
    for row, parents in parts:
        for parent in parents:
            if parent not in transcripts:
                raise ValueError('Unresolved model block Parent')
            model = transcripts[parent]
            if (row['seq'], row['strand']) != (model['seq'], model['strand']) or not model['start'] <= row['start'] < row['end'] <= model['end']:
                raise ValueError(f"Model block differs from transcript coordinates/strand: {parent}; "
                    f"{row['kind']} {row['seq']}:{row['start']+1}-{row['end']} ({row['strand']}); "
                    f"transcript {model['seq']}:{model['start']+1}-{model['end']} ({model['strand']})")
            model['exons' if row['kind'] == 'exon' else 'cds'].append((row['start'], row['end']))
    for tid, row in transcripts.items():
        if row['gene_id'] not in genes or not row['exons']:
            raise ValueError(f'{tid}: unresolved gene or missing transcript exons')
        gene = genes[row['gene_id']]
        if (row['seq'], row['strand']) != (gene['seq'], gene['strand']) or not gene['start'] <= row['start'] < row['end'] <= gene['end']:
            raise ValueError('Transcript differs from gene bounds/strand')
    if not transcripts:
        raise ValueError('No models to score')
    return transcripts

def rna_hits(path, lengths):
    transcripts, exons = {}, []
    for row in features(path, lengths, gtf=True):
        tid = row['attrs'].get('transcript_id')
        if row['kind'] == 'transcript':
            if not tid or tid in transcripts:
                raise ValueError('Duplicate/missing RNA transcript ID')
            transcripts[tid] = row | dict(identifier=tid, blocks=[])
        elif row['kind'] == 'exon':
            exons.append((tid, row))
    for tid, exon in exons:
        if tid not in transcripts:
            raise ValueError('RNA exon has unresolved transcript')
        hit = transcripts[tid]
        if (exon['seq'], exon['strand']) != (hit['seq'], hit['strand']) or not hit['start'] <= exon['start'] < exon['end'] <= hit['end']:
            raise ValueError('RNA exon outside transcript bounds/strand')
        hit['blocks'].append((exon['start'], exon['end']))
    for hit in transcripts.values():
        if not hit['blocks']:
            raise ValueError('RNA transcript has no exons')
        yield hit
