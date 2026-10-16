#!/usr/bin/env python3
"""Adapt and validate native InGenAnnot UTR-only refinement; no UTR prediction."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from evidence_aed import features, models, merged, overlap_length, rna_hits
from prepare_protein_evidence import checksum, records
from reference_repeats import fasta


def read(path): return json.loads(Path(path).read_text())
def dump(path,data): Path(path).write_text(json.dumps(data,indent=2)+'\n')


def subtract(blocks,excluded):
    result=[]
    for start,end in merged(blocks):
        parts=[(start,end)]
        for a,b in merged(excluded):
            parts=[piece for s,e in parts for piece in ((s,min(e,a)),(max(s,b),e)) if piece[1]>piece[0]]
        result.extend(parts)
    return merged(result)


def prepare(gff,validation,identity,rna,metadata,bams,indices,output):
    state=read(validation); reference=read(identity)
    if state['status']!='ready' or checksum(gff)!=state['gff3_sha256'] or state['reference_identity']!=reference['reference_identity']:
        raise ValueError('UTR inputs differ from accepted predictions/reference')
    output=Path(output); output.mkdir(exist_ok=True)
    rows=read(metadata)
    if len(rows)!=len(bams) or len(rows)!=len(indices): raise ValueError('UTR BAM metadata count differs')
    if rna and not rows: raise ValueError('RNA refinement requires BAMs for ranking')
    for row,bam,index in zip(rows,bams,indices):
        if row['layout'] not in ('paired','single') or row['strandedness'] not in ('reverse','unstranded'):
            raise ValueError('Pinned InGenAnnot ranking supports reverse/unstranded BAMs; forward/auto requires a validated adapter')
        if not Path(index).is_file() or not Path(bam).is_file(): raise ValueError('Missing BAM/index')
    lengths={c['id']:c['length'] for c in reference['contigs']}
    resolved=ambiguous=0
    if rna:
        for hit in rna_hits(rna,lengths):
            if hit['strand']=='.': ambiguous+=1
            else: resolved+=1
        if not resolved: raise ValueError('No resolved-strand RNA transcripts available for UTR refinement')
    dump(output/'inputs.json',dict(status='ready',model_sha256=checksum(gff),reference_identity=reference['reference_identity'],
        rna_sha256=checksum(rna) if rna else None,rna_transcripts=resolved,rna_ambiguous_strand_excluded=ambiguous,readsets=rows))


def finalize(original,native,identity,output):
    """Preserve BRAKER's explicit codons, excluding stop bases from native UTRs."""
    lengths={c['id']:c['length'] for c in read(identity)['contigs']}
    raw=list(features(native,lengths)); before=list(features(original,lengths))
    original_models=models(original,lengths)
    native_transcripts={r['attrs']['ID']:r for r in raw if r['kind'] in ('mRNA','transcript')}
    if set(native_transcripts)!=set(original_models):
        raise ValueError('Native UTR refinement changed transcript IDs')
    native_exons={tid:[] for tid in original_models}
    def coding_signature(rows):
        coding={tid:[] for tid in original_models}
        for row in rows:
            if row['kind']=='CDS':
                for tid in row['attrs']['Parent'].split(','):
                    if tid not in coding: raise ValueError('Native CDS has unresolved transcript parent')
                    coding[tid].append((row['seq'],row['strand'],row['start'],row['end'],row['phase']))
        return {tid:sorted(parts) for tid,parts in coding.items()}
    # Never allow restoring a rejected exon proposal to conceal a coding change.
    if coding_signature(before)!=coding_signature(raw):
        raise ValueError('Native UTR refinement changed CDS coordinates, strand or phases')
    rejected={}
    for row in raw:
        if row['kind'] in ('gene','mRNA','transcript'): continue
        for tid in row['attrs'].get('Parent','').split(','):
            if tid not in native_transcripts: raise ValueError('Native UTR feature has unresolved transcript parent')
            parent=native_transcripts[tid]; original_model=original_models[tid]
            if (row['seq'],row['strand'])!=(parent['seq'],parent['strand']) or (parent['seq'],parent['strand'])!=(original_model['seq'],original_model['strand']):
                rejected[tid]='Native UTR candidate has incompatible feature/transcript contig or strand; retained original BRAKER transcript'
            if row['kind']=='exon': native_exons[tid].append((row['start'],row['end']))
    for tid,old in original_models.items():
        exons=native_exons[tid]
        parent=native_transcripts[tid]
        if parent['attrs'].get('Parent')!=old['gene_id']:
            raise ValueError('Native UTR refinement changed gene/transcript relationships')
        reason=None
        if sorted(exons)!=sorted(old['exons']):
            if not parent['attrs'].get('utr_refine_evidence'):
                reason='Native exon change lacks UTR evidence provenance'
            elif overlap_length(exons,old['exons'])!=sum(e-s for s,e in merged(old['exons'])):
                reason='Native UTR candidate shortens existing exons'
            elif subtract(old['exons'],old['cds']):
                reason='Native UTR candidate changes existing UTRs despite onlynew'
            elif overlap_length(subtract(exons,old['exons']),[(min(s for s,e in old['cds']),max(e for s,e in old['cds']))]):
                reason='Native UTR candidate adds exonic bases inside the CDS span'
        if reason and tid not in rejected:
            rejected[tid]=reason+'; retained original BRAKER transcript'
    if rejected:
        restored=[]
        for row in raw:
            if row['kind']=='gene': restored.append(row); continue
            if row['kind'] in ('mRNA','transcript'):
                if row['attrs']['ID'] not in rejected: restored.append(row)
                continue
            parents=[p for p in row['attrs']['Parent'].split(',') if p not in rejected]
            if parents: restored.append(dict(row,attrs=dict(row['attrs'],Parent=','.join(parents))))
        for row in before:
            if row['kind'] in ('mRNA','transcript') and row['attrs']['ID'] in rejected:
                restored.append(row)
            elif row['kind'] in ('exon','CDS','five_prime_UTR','three_prime_UTR'):
                parents=[p for p in row['attrs']['Parent'].split(',') if p in rejected]
                if parents: restored.append(dict(row,attrs=dict(row['attrs'],Parent=','.join(parents))))
        raw=restored
    codons=[r for r in before if r['kind'] in ('start_codon','stop_codon')]
    stops={}
    for row in codons:
        if row['kind']=='stop_codon':
            for tid in row['attrs']['Parent'].split(','): stops.setdefault(tid,[]).append((row['start'],row['end']))
    adjusted=[]
    for row in raw:
        if row['kind'] in ('start_codon','stop_codon'): continue
        if row['kind'] in ('five_prime_UTR','three_prime_UTR'):
            tids=row['attrs']['Parent'].split(',')
            for tid in tids:
                for i,(s,e) in enumerate(subtract([(row['start'],row['end'])],stops.get(tid,[]))):
                    attrs=dict(row['attrs'],Parent=tid,ID=row['attrs']['ID']+f'.{i+1}')
                    adjusted.append(dict(row,start=s,end=e,attrs=attrs))
        else: adjusted.append(row)
    adjusted.extend(codons)
    transcripts={r['attrs']['ID']:r for r in adjusted if r['kind'] in ('mRNA','transcript')}
    genes={r['attrs']['ID']:r for r in adjusted if r['kind']=='gene'}
    for row in adjusted:
        if row['kind'] in ('gene','mRNA','transcript'): continue
        for tid in row['attrs']['Parent'].split(','):
            t=transcripts[tid]; t['start']=min(t['start'],row['start']); t['end']=max(t['end'],row['end'])
    gene_children={gid:[] for gid in genes}
    for t in transcripts.values(): gene_children[t['attrs']['Parent']].append(t)
    for gid,g in genes.items():
        children=gene_children[gid]
        if not children: raise ValueError('Native UTR output has a gene without transcripts')
        g['start']=min(t['start'] for t in children); g['end']=max(t['end'] for t in children)
    # Stable gene/transcript grouping is compatible with the existing filters.
    def gene_of(row):
        if row['kind']=='gene': return row['attrs']['ID']
        if row['kind'] in ('mRNA','transcript'): return row['attrs']['Parent']
        return transcripts[row['attrs']['Parent'].split(',')[0]]['attrs']['Parent']
    adjusted.sort(key=lambda r:(gene_of(r),0 if r['kind']=='gene' else 1 if r['kind'] in ('mRNA','transcript') else 2,r['start'],r['end'],r['kind']))
    with Path(output).open('w') as handle:
        handle.write('##gff-version 3\n')
        for r in adjusted:
            attrs=';'.join(k+'='+v for k,v in r['attrs'].items())
            handle.write(f"{r['seq']}\tingenannot\t{r['kind']}\t{r['start']+1}\t{r['end']}\t.\t{r['strand']}\t{r['phase']}\t{attrs}\n")
    dump(Path(output).with_suffix('.rejections.json'),rejected)


def sequences(path):
    result={}
    with Path(path).open() as handle:
        for h,s in records(handle):
            tid=h.split()[0]
            if tid in result: raise ValueError('Duplicate extracted sequence ID')
            result[tid]=s.upper()
    return result


def validate(original,refined,identity,genome,inputs,old_cds,new_cds,old_proteins,new_proteins,output,skipped=False,rejections=None):
    accepted=read(inputs); reference=read(identity)
    if checksum(original)!=accepted['model_sha256'] or reference['reference_identity']!=accepted['reference_identity']:
        raise ValueError('UTR annotation/reference changed since preparation')
    contigs=[dict(id=h.split()[0],length=len(s),uppercase_sha256=hashlib.sha256(s.upper().encode('ascii')).hexdigest()) for h,s in fasta(genome)]
    if contigs!=reference['contigs']: raise ValueError('UTR genome differs from accepted reference')
    lengths={c['id']:c['length'] for c in reference['contigs']}
    before,after=models(original,lengths),models(refined,lengths)
    rejected=read(rejections) if rejections else {}
    if not set(rejected)<=set(before): raise ValueError('UTR rejection refers to unknown transcript')
    if set(before)!=set(after): raise ValueError('UTR refinement changed transcript IDs')
    def coding(path):
        result={tid:[] for tid in before}
        for row in features(path,lengths):
            if row['kind']=='CDS':
                for tid in row['attrs']['Parent'].split(','): result[tid].append((row['start'],row['end'],row['phase']))
        return {k:sorted(v) for k,v in result.items()}
    if coding(original)!=coding(refined): raise ValueError('UTR refinement changed CDS coordinates or phases')
    old_rows=list(features(original,lengths)); new_rows=list(features(refined,lengths))
    if {r['attrs']['ID'] for r in old_rows if r['kind']=='gene'}!={r['attrs']['ID'] for r in new_rows if r['kind']=='gene'}:
        raise ValueError('UTR refinement changed gene IDs')
    ids={}; stops={}; utrs={tid:{'five_prime_UTR':[],'three_prime_UTR':[]} for tid in before}
    for row in old_rows:
        if row['kind']=='stop_codon':
            for tid in row['attrs']['Parent'].split(','): stops.setdefault(tid,[]).append((row['start'],row['end']))
    for row in new_rows:
        ident=row['attrs'].get('ID')
        if not ident: raise ValueError('Refined feature missing ID')
        if ident in ids and not (row['kind']=='CDS' and ids[ident]=='CDS'): raise ValueError('Duplicate refined feature ID')
        ids[ident]=row['kind']
    for row in new_rows:
        parents=row['attrs'].get('Parent','').split(',')
        if row['kind']!='gene' and (not parents[0] or any(p not in ids for p in parents)): raise ValueError('Unresolved refined Parent')
        if row['kind'] not in ('gene','mRNA','transcript'):
            for tid in parents:
                if tid not in after: raise ValueError('Refined feature parent is not a transcript')
                parent=after[tid]
                if (row['seq'],row['strand'])!=(parent['seq'],parent['strand']) or not parent['start']<=row['start']<row['end']<=parent['end']:
                    raise ValueError('Refined feature outside parent bounds/strand')
        if row['kind'] in ('five_prime_UTR','three_prime_UTR'):
            for tid in parents: utrs[tid][row['kind']].append((row['start'],row['end']))
    summary=[]
    for tid,a in before.items():
        b=after[tid]
        if tid in rejected and (sorted(a['exons'])!=sorted(b['exons']) or a['attrs']!=b['attrs']):
            raise ValueError('Rejected UTR candidate did not retain original transcript')
        if any(a[k]!=b[k] for k in ('gene_id','seq','strand')): raise ValueError('UTR refinement changed gene/reference/strand')
        if overlap_length(b['exons'],a['exons'])!=sum(e-s for s,e in merged(a['exons'])): raise ValueError('UTR refinement shortened existing exons')
        if overlap_length(b['exons'],a['cds'])!=sum(e-s for s,e in merged(a['cds'])): raise ValueError('Refined CDS outside exons')
        if overlap_length(subtract(b['exons'],a['exons']),[(min(s for s,e in a['cds']),max(e for s,e in a['cds']))]):
            raise ValueError('UTR refinement added exonic bases inside CDS span')
        if subtract(a['exons'],a['cds']) and sorted(a['exons'])!=sorted(b['exons']): raise ValueError('Existing UTRs changed despite onlynew')
        if sorted(a['exons'])!=sorted(b['exons']) and not b['attrs'].get('utr_refine_evidence'):
            raise ValueError(f'UTR extension lacks native evidence provenance: {tid}')
        lower,upper=min(s for s,e in a['cds']),max(e for s,e in a['cds'])
        low=[(s,min(e,lower)) for s,e in b['exons'] if s<lower]
        high=[(max(s,upper),e) for s,e in b['exons'] if e>upper]
        expected={'five_prime_UTR':low if a['strand']=='+' else high,'three_prime_UTR':high if a['strand']=='+' else low}
        for kind,blocks in expected.items():
            if not skipped and tid not in rejected and merged(utrs[tid][kind])!=subtract(blocks,stops.get(tid,[])): raise ValueError('UTR orientation/coordinates differ from exons or include stop codon')
        summary.append(dict(transcript_id=tid,gene_id=a['gene_id'],changed=sorted(a['exons'])!=sorted(b['exons']),
            rejected_candidate_reason=rejected.get(tid,'.'),
            evidence_id=b['attrs'].get('utr_refine_evidence','.'),
            five_prime_utr_bp=sum(e-s for s,e in subtract(expected['five_prime_UTR'],stops.get(tid,[]))),
            three_prime_utr_bp=sum(e-s for s,e in subtract(expected['three_prime_UTR'],stops.get(tid,[])))))
    if sequences(old_cds)!=sequences(new_cds) or sequences(old_proteins)!=sequences(new_proteins): raise ValueError('UTR refinement changed extracted CDS/protein sequences')
    if set(sequences(new_cds))!=set(before) or set(sequences(new_proteins))!=set(before): raise ValueError('Extracted sequence IDs differ from model IDs')
    output=Path(output); output.mkdir(exist_ok=True)
    with (output/'utr_changes.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(summary[0]),delimiter='\t'); w.writeheader(); w.writerows(summary)
    status=dict(status='ready',stage='utr_refinement',scope='Native ranked RNA UTR refinement; no AED filtering',
        skipped=skipped,transcripts=len(before),changed_transcripts=sum(r['changed'] for r in summary),
        cds_coordinates_and_phases_unchanged=True,cds_sequences_unchanged=True,protein_sequences_unchanged=True,
        baseline_protein_qc_reusable=True,original_model_sha256=checksum(original),gff3_sha256=checksum(refined),
        reference_identity=reference['reference_identity'],readsets=accepted['readsets'])
    status['rna_transcripts']=accepted.get('rna_transcripts',0)
    status['rejected_utr_candidates']=rejected
    status['rna_ambiguous_strand_excluded']=accepted.get('rna_ambiguous_strand_excluded',0)
    dump(output/'utr_status.json',status)
    (output/'README.md').write_text('# UTR refinement\n\nStatus: **ready**\n\n'+
        ('Refinement skipped: disabled or no RNA evidence.\n\n' if skipped else 'Native clusterize, BAM-based isoform ranking and utr_refine --utr_mode rank --onlynew.\n\n')+
        f"Changed transcripts: {status['changed_transcripts']}. CDS coordinates, phases and extracted CDS/protein sequences are unchanged. "
        f"Rejected UTR candidates: {len(rejected)}; original BRAKER transcripts retained. "
        'Original BRAKER models remain preserved. Separately annotated stop codons are retained and excluded from UTR lengths.\n\n'
        '[Changes](utr_changes.tsv) | [Validation](utr_status.json)\n\n'
        'Review extensions before promotion. Protein BUSCO/OMArk inputs are unchanged; final AGAT/AED assessment follows. No models were filtered.\n')
    return status


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare')
    for arg in ('gff','validation','identity','metadata','outdir'): a.add_argument('--'+arg,required=True)
    a.add_argument('--rna'); a.add_argument('--bams',nargs='*',default=[]); a.add_argument('--indices',nargs='*',default=[])
    a=sub.add_parser('finalize')
    for arg in ('original','native','identity','output'): a.add_argument('--'+arg,required=True)
    a=sub.add_parser('validate')
    for arg in ('original','refined','identity','genome','inputs','old-cds','new-cds','old-proteins','new-proteins','outdir'): a.add_argument('--'+arg,required=True)
    a.add_argument('--skipped',action='store_true')
    a.add_argument('--rejections')
    a=p.parse_args()
    if a.action=='prepare': prepare(a.gff,a.validation,a.identity,a.rna,a.metadata,a.bams,a.indices,a.outdir)
    elif a.action=='finalize': finalize(a.original,a.native,a.identity,a.output)
    else: validate(a.original,a.refined,a.identity,a.genome,a.inputs,a.old_cds,a.new_cds,a.old_proteins,a.new_proteins,a.outdir,a.skipped,a.rejections)


if __name__=='__main__': main()
