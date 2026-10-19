#!/usr/bin/env python3
"""Join structural metrics, repeat flags and evidence QC without model selection."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

from prepare_protein_evidence import checksum


def write(path,rows):
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({k:'.' if v is None else v for k,v in row.items()})


def load_table(path,key):
    with Path(path).open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    if not rows or any(not r.get(key) for r in rows) or len({r[key] for r in rows})!=len(rows):
        raise ValueError('Missing/duplicate assessment evidence IDs')
    return {r[key]:r for r in rows}


def evidence(status_path,table_path,metrics,key,expected):
    if not status_path and not table_path:return None,{}
    if not status_path or not table_path:raise ValueError('Evidence table requires its validation status')
    state=json.loads(Path(status_path).read_text())
    if (state['status']!='ready' or state['model_sha256']!=metrics['model_sha256']
        or state['reference_identity']!=metrics['reference_identity']):raise ValueError('Evidence/model/reference mismatch')
    rows=load_table(table_path,key)
    if set(rows)!=set(expected):raise ValueError('Evidence IDs differ from accepted annotation')
    for k,row in rows.items():
        if key=='transcript_id' and row.get('gene_id')!=expected[k]['gene_id']:raise ValueError('Evidence transcript gene parent mismatch')
    return state,rows


def report(aid,metrics,outdir,aed_status=None,aed_scores=None,junction_status=None,junction_scores=None,te_status=None,te_scores=None):
    metrics=json.loads(Path(metrics).read_text())
    if metrics['status']!='ready' or metrics['assembly_id']!=aid:raise ValueError('Model metrics not ready for assembly')
    tx={r['transcript_id']:r for r in metrics['metrics']};by_gene=defaultdict(list)
    if len(tx)!=len(metrics['metrics']):raise ValueError('Duplicate model metric transcript')
    for tid,row in tx.items():by_gene[row['gene_id']].append(row)
    astate,aed=evidence(aed_status,aed_scores,metrics,'transcript_id',tx)
    jstate,junctions=evidence(junction_status,junction_scores,metrics,'transcript_id',tx)
    tstate,te=evidence(te_status,te_scores,metrics,'gene_id',by_gene)
    if astate and astate.get('stage')!='native_aed':raise ValueError('Expected production native AED')
    if jstate and jstate.get('stage')!='junction_qc':raise ValueError('Expected junction QC')
    if tstate:
        if tstate.get('stage')!='te_flags' or checksum(te_scores)!=tstate['report_sha256']:raise ValueError('TE report changed or wrong stage')
    rows=[];gene_rows=[]
    for tid,metric in sorted(tx.items()):
        row=dict(metric);gid=row['gene_id']
        row.update(aed_assessed=bool(astate),junction_qc_enabled=bool(jstate),te_assessed=bool(tstate))
        for key in ('rna_aed','protein_aed','rna_evidence_id','protein_evidence_id','rna_penalty'):
            row[key]=aed.get(tid,{}).get(key,'.')
        if astate:
            for kind in ('rna','protein'):
                value=row[kind+'_aed'];available=astate['evidence'][kind]['available']
                if (value!='.')!=available:raise ValueError('AED availability mismatch')
                if available and (not math.isfinite(float(value)) or not 0<=float(value)<=1):raise ValueError('Invalid AED value')
        jrow=junctions.get(tid)
        if jrow and int(jrow['all_junctions'])!=row['n_introns']:raise ValueError('Junction/model intron count mismatch')
        row['junction_status']=jrow['status'] if jrow else 'disabled'
        for prefix in ('all','coding','utr'):
            for suffix in ('junctions','pass','fail','not_observed','unresolved_strand','pass_fraction'):
                key=prefix+'_'+suffix;row['junction_'+key]=jrow[key] if jrow else None
        trow=te.get(gid)
        row['te_gene_flag']=trow['te_flagged'] if trow else None
        row['te_evidence']=trow['te_evidence'] if trow else None
        for kind in ('gene','exon','cds'):row['te_'+kind+'_repeat_fraction']=float(trow[kind+'_repeat_fraction']) if trow else None
        hit=tstate['transcript_hits'].get(tid) if tstate else None
        row['te_transcript_homology']=bool(hit) if tstate and tstate['search']['assessed'] else None
        for key in ('hit','pident','qcov','evalue'):row['te_transcript_'+key]=hit[key] if hit else None
        rows.append(row)
    ranks={'PASS':0,'WARN':1,'FAIL':2}
    for gid,group in sorted(by_gene.items()):
        flags=[r['structure_flag'] for r in group]
        gene_rows.append(dict(gene_id=gid,transcripts=len(group),worst_structure_flag=max(flags,key=ranks.get),
            best_structure_flag=min(flags,key=ranks.get),n_pass=flags.count('PASS'),n_warn=flags.count('WARN'),n_fail=flags.count('FAIL'),
            min_protein_length_aa=min(r['protein_length_aa'] for r in group),max_protein_length_aa=max(r['protein_length_aa'] for r in group),
            te_flagged=te[gid]['te_flagged'] if te else None,te_evidence=te[gid]['te_evidence'] if te else None,
            reasons=';'.join(r['transcript_id']+':'+r['structure_reasons'] for r in group if r['structure_reasons']!='.') or '.'))
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    write(out/'transcript_assessment.tsv',rows);write(out/'gene_assessment.tsv',gene_rows)
    status=dict(status='ready',stage='model_assessment',assembly_id=aid,model_sha256=metrics['model_sha256'],reference_identity=metrics['reference_identity'],
        transcripts=len(rows),genes=len(gene_rows),models_removed=0,accepted_models_modified=False,
        structure_flags=dict(Counter(r['structure_flag'] for r in rows)),aed_enabled=bool(astate),junction_qc_enabled=bool(jstate),te_enabled=bool(tstate),
        te_flagged_genes=tstate['flagged_genes'] if tstate else None,
        source_metrics='Structural metrics and severity rules adapted from MERGE_FUNCTIONAL_ANNOTATIONS; gene best/worst aggregation from STRUCTURAL_QC; functional support not yet assessed',
        scope='Diagnostic flags on every accepted transcript/gene; no filtering, renaming or combined confidence score',
        sources=dict(aed=checksum(aed_scores) if aed_scores else None,junctions=checksum(junction_scores) if junction_scores else None,te=checksum(te_scores) if te_scores else None),
        outputs={name:checksum(out/name) for name in ('transcript_assessment.tsv','gene_assessment.tsv')})
    (out/'assessment_status.json').write_text(json.dumps(status,indent=2)+'\n')
    (out/'README.md').write_text('# Model assessment before selection\n\nStatus: **ready**\n\n'
        f"Genes: {len(gene_rows)}. Coding transcripts: {len(rows)}. Models removed: **0**.\n\n"
        '[Transcript assessment](transcript_assessment.tsv) | [Gene summary](gene_assessment.tsv) | [Provenance](assessment_status.json)\n\n'
        'Structure flags retain the original sequence/exon/splice checks. FAIL/WARN describe diagnostic severity, not pipeline failure or exclusion. '
        'Functional support and OMArk-dependent quality rules are deferred until their inputs are connected. '
        'Absence of a terminal protein * or an explicit codon feature does not prove an incomplete gene.\n\n'
        'Repeat overlap is the union of masked intervals over gene span; CDS/exonic union fractions are also shown. '
        'Homology to the pre-filter repeat library is repeat-associated evidence, not definitive TE identity; that library can include host-gene sequences. '
        'A gene flag can come from another isoform: transcript homology is shown separately. '
        'Missing/disabled evidence is unavailable, never replaced with a poor score. No combined AED/confidence score is created.\n\n'
        f'[Original accepted models](../../final/{aid}/README.md) | [TE annotations](te/{aid}.te_flagged.gff3) (when enabled)\n')
    return status


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','metrics','outdir'):p.add_argument('--'+key,required=True)
    for key in ('aed-status','aed-scores','junction-status','junction-scores','te-status','te-scores'):p.add_argument('--'+key)
    a=p.parse_args();report(a.assembly_id,a.metrics,a.outdir,a.aed_status,a.aed_scores,a.junction_status,a.junction_scores,a.te_status,a.te_scores)
