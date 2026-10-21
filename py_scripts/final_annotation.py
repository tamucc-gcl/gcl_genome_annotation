#!/usr/bin/env python3
"""Publish validated final models and report native evidence-specific AED."""
import argparse
import csv
import json
from pathlib import Path
import statistics

from ingenannot_aed import assess
from prepare_protein_evidence import checksum


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')



def aed(gff,annotated,identity,inputs,cds_only,penalty,outdir):
    before,scores,normalized,accepted=assess(gff,annotated,identity,inputs)
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for tid,model in sorted(before.items()):
        score=scores[tid]
        rows.append(dict(gene_id=model['gene_id'],transcript_id=tid,rna_aed=score['rna'],protein_aed=score['protein'],
            rna_evidence_id=score['rna_evidence_id'] if score['rna'] is not None else None,
            protein_evidence_id=score['protein_evidence_id'] if score['protein'] is not None else None,
            rna_penalty=score['rna_penalty'] if score['rna'] is not None else None))
    with (outdir/'aed_scores.tsv').open('w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({k:'.' if v is None else v for k,v in row.items()})
    status=dict(status='ready',stage='native_aed',model_sha256=checksum(gff),reference_identity=accepted['reference_identity'],
        evidence=accepted['evidence'],transcripts=len(rows),genes=len({r['gene_id'] for r in rows}),
        profile=dict(tool='InGenAnnot',version='0.0.20',rna_cds_only=cds_only,rna_overflow_penalty=penalty,protein_cds_only=True),
        models_unchanged=True,transcript_bounds_normalized_to_exons=normalized,
        native_exon_record_normalizations=accepted['native_exon_record_normalizations'],
        scope='Native RNA and protein AED on the accepted final annotation; no combined score or model filtering')
    status['metrics']={}
    for kind in ('rna','protein'):
        values=[r[kind+'_aed'] for r in rows if r[kind+'_aed'] is not None]
        status['metrics'][kind]=dict(available=len(values),median=statistics.median(values) if values else None,
            mean=statistics.mean(values) if values else None,score_equals_zero=sum(v==0 for v in values),
            score_equals_one=sum(v==1 for v in values))
    dump(outdir/'aed_status.json',status)
    (outdir/'README.md').write_text('# Native annotation evidence assessment\n\nStatus: **ready**\n\n'
        f"Transcripts: {len(rows)}. Separate RNA and protein AED; missing evidence is `.`, not AED 1.\n\n"
        f"RNA CDS-only profile: {cds_only}. RNA overflow penalty: {penalty}. Protein assessment uses CDS. "
        'CDS-only RNA assessment measures coding support, not UTR boundary quality. '
        'Lower AED means closer agreement with selected evidence; it is not a probability of correctness.\n\n'
        '[Scores](aed_scores.tsv) | [Settings and validation](aed_status.json) | [Native output](native/native.aed.gff3)\n\n'
        '[RNA junction support](../../../qc/junctions/'+str(accepted.get('assembly_id') or Path(gff).stem)+'/README.md) (when enabled).\n\n'
        f"Repeated identical exon records differ for {len(accepted['native_exon_record_normalizations'])} transcripts in the native export; counts are audited in the status JSON. "
        'Exact unique exon boundaries, CDS records and coding phases are preserved.\n\n'
        'Assessment does not replace or edit the accepted final GFF3. Filtering is disabled by default; '
        'source-specific filtering is available in the downstream selection stage.\n')
    return status


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('aed')
    for key in ('gff','annotated','identity','inputs','outdir'):p.add_argument('--'+key,required=True)
    p.add_argument('--cds-only',choices=['true','false'],required=True);p.add_argument('--penalty',type=float,required=True)
    a=parser.parse_args()
    aed(a.gff,a.annotated,a.identity,a.inputs,a.cds_only=='true',a.penalty,a.outdir)
