#!/usr/bin/env python3
"""Structural metrics moved upstream from MERGE_FUNCTIONAL_ANNOTATIONS.

Retains its sequence/exon/splice checks without requiring functional searches.
Uses validated parents, unique blocks and correct internal-stop handling.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from urllib.parse import quote

from evidence_aed import models, features, merged
from prepare_protein_evidence import records, checksum
from reference_repeats import fasta


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def sequences(path):
    result={}
    with Path(path).open() as handle:
        for header,seq in records(handle):
            tid=header.split()[0]
            if tid in result or not seq:raise ValueError('Duplicate/empty sequence ID')
            result[tid]=seq.upper()
    return result


def revcomp(seq):return seq.translate(str.maketrans('ACGTNRYSWKMBDHV','TGCANYRSWMKVHDB'))[::-1]


def splice_class(donor,acceptor):
    # Same GT-AG, GC-AG and AT-AC categories as the original functional merger.
    return donor+acceptor if (donor,acceptor) in [('GT','AG'),('GC','AG'),('AT','AC')] else 'non_canonical'


def prepare(aid,gff,proteins,cds,genome,identity,validation,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    state=json.loads(Path(validation).read_text());ref=json.loads(Path(identity).read_text())
    if (state['status']!='ready' or state.get('stage')!='final_annotation' or state.get('assembly_id')!=aid
        or state['reference_identity']!=ref['reference_identity']):raise ValueError('Assessment requires accepted final models/reference')
    for path,key in [(gff,'gff3_sha256'),(proteins,'proteins_sha256'),(cds,'cds_sha256')]:
        if checksum(path)!=state[key]:raise ValueError('Accepted model/sequence file changed')
    genome_seqs={h.split()[0]:seq.upper() for h,seq in fasta(genome)}
    observed=[dict(id=k,length=len(v),uppercase_sha256=hashlib.sha256(v.encode()).hexdigest()) for k,v in genome_seqs.items()]
    if observed!=ref['contigs']:raise ValueError('Assessment genome differs from accepted reference')
    lengths={c['id']:c['length'] for c in ref['contigs']};tx=models(gff,lengths)
    protein=sequences(proteins);coding=sequences(cds)
    if set(tx)!=set(protein) or set(tx)!=set(coding):raise ValueError('GFF/protein/CDS transcript IDs differ')
    gene_rows={r['attrs']['ID']:r for r in features(gff,lengths) if r['kind']=='gene'}
    codons=defaultdict(set)
    for row in features(gff,lengths):
        if row['kind'] in ('start_codon','stop_codon'):
            for tid in row['attrs'].get('Parent','').split(','):
                if tid not in tx:raise ValueError('Unresolved codon Parent')
                codons[tid].add(row['kind'])
    metric_rows=[];gene_blocks=defaultdict(lambda:defaultdict(list))
    for tid,m in sorted(tx.items()):
        exons=sorted(set(m['exons']));blocks=sorted(set(m['cds']))
        if not blocks:raise ValueError('Coding model has no CDS')
        if any(l[1]>=r[0] for l,r in zip(exons,exons[1:])):raise ValueError('Overlapping/adjacent exon blocks')
        if sum(e-s for s,e in blocks)!=sum(e-s for s,e in merged(blocks)):raise ValueError('Overlapping CDS blocks')
        if any(not any(x<=s<e<=y for x,y in exons) for s,e in blocks):raise ValueError('CDS block outside exons')
        introns=[(l[1],r[0]) for l,r in zip(exons,exons[1:])];classes=Counter()
        for s,e in introns:
            seq=genome_seqs[m['seq']]
            if e-s<2:classes['non_canonical']+=1;continue
            donor,acceptor=(seq[s:s+2],seq[e-2:e]) if m['strand']=='+' else (revcomp(seq[e-2:e]),revcomp(seq[s:s+2]))
            classes[splice_class(donor,acceptor)]+=1
        seq=protein[tid];body=seq[:-1] if seq.endswith('*') else seq
        n_noncan=classes['non_canonical'];n_can=len(introns)-n_noncan
        row=dict(gene_id=m['gene_id'],transcript_id=tid,contig=m['seq'],strand=m['strand'],
            gene_span_bp=gene_rows[m['gene_id']]['end']-gene_rows[m['gene_id']]['start'],transcript_span_bp=m['end']-m['start'],
            protein_length_aa=len(body),has_start_M='Y' if seq.startswith('M') else 'N',
            has_internal_stop='Y' if '*' in body else 'N',has_terminal_stop_symbol='Y' if seq.endswith('*') else 'N',
            has_explicit_start_codon='Y' if 'start_codon' in codons[tid] else 'N',
            has_explicit_stop_codon='Y' if 'stop_codon' in codons[tid] else 'N',
            cds_sequence_bp=len(coding[tid]),cds_length_bp=sum(e-s for s,e in blocks),n_exons=len(exons),
            is_single_exon='Y' if len(exons)==1 else 'N',min_exon_len=min(e-s for s,e in exons),
            n_introns=len(introns),min_intron_len=min((e-s for s,e in introns),default=None),
            n_short_introns_lt60=sum(e-s<60 for s,e in introns),n_canonical_introns=n_can,n_noncanonical_introns=n_noncan,
            pct_canonical_splice=100*n_can/len(introns) if introns else None,
            donor_acceptor_summary=','.join(f'{k}:{v}' for k,v in sorted(classes.items())) or '.')
        # Original structural severity rules, with complete reasons retained.
        severe=[];warnings=[]
        if row['has_internal_stop']=='Y':severe.append('internal_stop_codon')
        if row['protein_length_aa']<30:severe.append('very_short_protein')
        if len(introns)>=2 and n_noncan/len(introns)>.5:severe.append('majority_noncanonical_splice')
        if row['has_start_M']=='N':warnings.append('no_start_met')
        if row['protein_length_aa']<100:warnings.append('short_protein')
        if row['n_short_introns_lt60']:warnings.append('short_introns_lt60')
        if introns and n_noncan/len(introns)>.1:warnings.append('low_canonical_splice_pct')
        row.update(structure_flag='FAIL' if severe else 'WARN' if warnings else 'PASS',structure_reasons=';'.join(severe+warnings) or '.')
        metric_rows.append(row)
        gene_blocks[m['gene_id']]['exon'].extend(exons);gene_blocks[m['gene_id']]['cds'].extend(blocks)
    with (out/'regions.bed').open('w') as handle:
        for gid,gene in sorted(gene_rows.items()):
            blocks=dict(gene=[(gene['start'],gene['end'])],exon=merged(gene_blocks[gid]['exon']),cds=merged(gene_blocks[gid]['cds']))
            for kind,parts in blocks.items():
                for s,e in parts:handle.write(f"{gene['seq']}\t{s}\t{e}\t{kind}|{quote(gid,safe='')}\t0\t{gene['strand']}\n")
    result=dict(status='ready',assembly_id=aid,stage='model_metrics',model_sha256=checksum(gff),reference_identity=ref['reference_identity'],
        proteins_sha256=checksum(proteins),cds_sha256=checksum(cds),genes=len(gene_rows),transcripts=len(tx),
        metrics=metric_rows,scope='Legacy structural/sequence/splice checks upstream of functional searches; diagnostic flags only')
    dump(out/'metrics.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','gff','proteins','cds','genome','identity','validation','outdir'):p.add_argument('--'+key,required=True)
    a=p.parse_args();prepare(a.assembly_id,a.gff,a.proteins,a.cds,a.genome,a.identity,a.validation,a.outdir)
