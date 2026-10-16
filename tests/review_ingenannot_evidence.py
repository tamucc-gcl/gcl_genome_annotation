#!/usr/bin/env python3
"""Extract a small coordinate-level review from completed AED outputs."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from evidence_aed import models, merged, nucleotide_aed, overlap_length, rna_hits
from prepare_protein_evidence import checksum


def table(path):
    with Path(path).open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    result={r['transcript_id']:r for r in rows}
    if len(result)!=len(rows): raise ValueError('Duplicate transcript IDs in score table')
    return result


def choose_cases(benchmark,baseline):
    selected={}
    def add(category,rows):
        genes=set(); count=0
        for row in rows:
            tid=row['transcript_id']
            if tid in selected or row['gene_id'] in genes: continue
            selected[tid]=category; genes.add(row['gene_id']); count+=1
            if count==2: break
    rows=list(benchmark.values())
    valid=[r for r in rows if r['ingenannot_rna_aed']!='.' and r['ingenannot_protein_aed']!='.']
    add('changed_rna_selection',sorted(
        (r for r in valid if r['rna_penalty']=='no' and float(r['rna_delta'])>.01
         and r['rna_evidence_id']!=baseline[r['transcript_id']]['rna_best_alignment_id']),
        key=lambda r:(-float(r['rna_delta']),r['transcript_id'])))
    add('rna_penalty',sorted((r for r in valid if r['rna_penalty']=='yes'),
        key=lambda r:(float(r['ingenannot_rna_aed']),r['transcript_id'])))
    add('weak_both_with_selected_evidence',sorted(
        (r for r in valid if .7<float(r['ingenannot_rna_aed'])<1 and .7<float(r['ingenannot_protein_aed'])<1),
        key=lambda r:(-(float(r['ingenannot_rna_aed'])+float(r['ingenannot_protein_aed'])),r['transcript_id'])))
    add('score_one_despite_previous_overlap',sorted(
        (r for r in valid if float(r['ingenannot_rna_aed'])==1 and float(r['current_rna_aed'])<1),
        key=lambda r:(float(r['current_rna_aed']),r['transcript_id'])))
    add('strong_agreement_control',sorted(
        (r for r in valid if float(r['ingenannot_rna_aed'])<.1 and float(r['ingenannot_protein_aed'])==0),
        key=lambda r:(float(r['ingenannot_rna_aed']),r['transcript_id'])))
    return selected


def compare_rna(model,hit,penalty):
    cds=merged(model['cds']); exons=merged(hit['blocks'])
    if not cds: raise ValueError('RNA CDS-only review requires coding models')
    common=overlap_length(cds,exons)
    cds_outside=sum(e-s for s,e in cds)-common
    inside=sum(max(0,min(e,cds[-1][1])-max(s,cds[0][0])) for s,e in exons)
    extra_inside=inside-common
    raw=nucleotide_aed(cds,exons)['aed']
    applied=raw<1 and (cds_outside>0 or extra_inside>0)
    model_introns={(a[1],b[0]) for a,b in zip(merged(model['exons']),merged(model['exons'])[1:])}
    rna_introns={(a[1],b[0]) for a,b in zip(exons,exons[1:])}
    return dict(evidence_id=hit['identifier'],seq=hit['seq'],strand=hit['strand'],exons=exons,
        raw_cds_aed=raw,penalty_applied=applied,penalized_aed=raw+(penalty if applied else 0),
        cds_bases_outside_rna=cds_outside,rna_extra_exonic_bases_inside_cds_span=extra_inside,
        model_introns=len(model_introns),rna_introns=len(rna_introns),
        matching_intron_boundaries=len(model_introns & rna_introns))


def selected_proteins(path,wanted,include_attributes=False):
    """Read grouped miniprot output once; keep only the requested chains."""
    hits={}; current=None
    if not wanted: return hits
    pattern=re.compile(r'(?:^|;)ID=([^;\s]+)')
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=9: continue
            if fields[2]=='mRNA':
                if wanted.issubset(hits): break
                match=pattern.search(fields[8]); identifier=match.group(1) if match else None
                current=identifier if identifier in wanted else None
                if current:
                    hits[current]=dict(evidence_id=current,seq=fields[0],strand=fields[6],blocks=[])
                    if include_attributes:
                        hits[current]['attributes']=dict(x.split('=',1) for x in fields[8].split(';') if '=' in x)
            elif fields[2]=='CDS' and current:
                if 'Parent='+current not in fields[8].split(';'):
                    raise ValueError('Unexpected miniprot CDS parent')
                hits[current]['blocks'].append((int(fields[3])-1,int(fields[4])))
    if set(hits)!=wanted or any(not h['blocks'] for h in hits.values()):
        raise ValueError('Requested protein chains missing from alignment source')
    return hits


def review(root,aid,output):
    root=Path(root); structural=root/'structural_annotation'
    folder=structural/'ingenannot_benchmark'/aid
    status=json.loads((folder/'benchmark_status.json').read_text())
    if status['status']!='ready' or not status['profile']['rna_cds_only']:
        raise ValueError('Review requires a completed CDS-only RNA benchmark')
    scores=table(folder/'benchmark_scores.tsv')
    baseline=table(structural/'evidence_aed'/aid/'aed_scores.tsv')
    if set(scores)!=set(baseline): raise ValueError('Baseline transcript IDs differ')
    for tid,row in scores.items():
        if any(row['current_'+k+'_aed']!=baseline[tid][k+'_aed'] for k in ('rna','protein')):
            raise ValueError('Baseline scores changed since benchmark')
    gff=structural/'braker'/aid/(aid+'.gff3')
    if checksum(gff)!=status['model_sha256']: raise ValueError('Annotation changed since benchmark')
    identity=json.loads((root/'filtered_assemblies'/aid/(aid+'.contigs.json')).read_text())
    if identity['reference_identity']!=status['reference_identity']: raise ValueError('Reference identity changed')
    lengths={c['id']:c['length'] for c in identity['contigs']}
    allmodels=models(gff,lengths)
    chosen=choose_cases(scores,baseline)
    if not chosen: raise ValueError('No review cases available')
    rna=structural/'braker'/aid/'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
    if checksum(rna)!=status['evidence']['rna']['source_sha256']: raise ValueError('RNA evidence changed')
    print('Reading retained RNA structures...',flush=True)
    hits=list(rna_hits(rna,lengths)); lookup={h['identifier']:h for h in hits}
    wanted={scores[t]['protein_evidence_id'] for t in chosen}-{'.','None',''}
    print('Extracting selected protein alignments (one streaming scan; no alignment)...',flush=True)
    proteins=selected_proteins(structural/'evidence_aed'/aid/'protein_alignment'/(aid+'.miniprot.gff'),wanted)
    cases=[]; penalty=status['profile']['rna_overflow_penalty']
    for tid,category in chosen.items():
        model=allmodels[tid]; row=scores[tid]; old=baseline[tid]
        candidates=[]
        for hit in hits:
            if hit['seq']==model['seq'] and hit['strand']==model['strand'] and hit['start']<model['end'] and hit['end']>model['start']:
                candidates.append(compare_rna(model,hit,penalty))
        best=min([1.0]+[h['penalized_aed'] for h in candidates])
        actual=float(row['ingenannot_rna_aed'])
        native=lookup.get(row['rna_evidence_id']); previous=lookup.get(old['rna_best_alignment_id'])
        protein=proteins.get(row['protein_evidence_id'])
        if protein and (protein['seq'],protein['strand'])!=(model['seq'],model['strand']):
            raise ValueError('Selected protein differs from model contig/strand')
        if protein: protein=dict(protein,raw_cds_aed=nucleotide_aed(model['cds'],protein['blocks'])['aed'])
        cases.append(dict(category=category,transcript_id=tid,gene_id=model['gene_id'],
            model=dict(seq=model['seq'],strand=model['strand'],exons=model['exons'],cds=model['cds']),
            scores=row,previous_rna=compare_rna(model,previous,penalty) if previous else None,
            native_rna=compare_rna(model,native,penalty) if native else None,
            top_rna_candidates=sorted(candidates,key=lambda h:(h['penalized_aed'],h['evidence_id']))[:3],
            recomputed_native_rna_aed=best,rna_score_matches=abs(best-actual)<=.000051,
            selected_protein=protein,protein_score_matches=abs(protein['raw_cds_aed']-float(row['ingenannot_protein_aed']))<=.000051 if protein else None))
    result=dict(stage='ingenannot_coordinate_review',assembly_id=aid,model_sha256=status['model_sha256'],
        coordinates='0-based half-open; introns compare assembled RNA exon boundaries, not independent read junction support',
        scope='Selected-case structural review; not a genome-wide accuracy validation or independent splice-read assessment',
        status='passed' if all(c['rna_score_matches'] and c['protein_score_matches'] is not False for c in cases) else 'disagreement',cases=cases)
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    (output/'review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"{result['status']}: {len(cases)} cases; {output/'review.json'}",flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs-root',required=True); p.add_argument('--assembly-id',required=True); p.add_argument('--outdir',required=True)
    a=p.parse_args(); review(a.inputs_root,a.assembly_id,a.outdir)
