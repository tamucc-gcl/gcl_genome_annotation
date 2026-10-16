#!/usr/bin/env python3
"""Validate and publish final models independently of AED reporting."""
import argparse
import json
from pathlib import Path
import shutil
from prepare_protein_evidence import checksum


def load(path):return json.loads(Path(path).read_text())
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def finalize(aid,prediction,validation,gff,utr_status,outdir,mode=None):
    original,utr=load(validation),load(utr_status)
    prediction,outdir=Path(prediction),Path(outdir)
    mode=mode or original.get('mode')
    if mode not in ('rna_protein','rna','protein','ab_initio') or original.get('mode',mode)!=mode:
        raise ValueError('Final annotation evidence mode differs from prediction')
    if original.get('assembly_id')!=aid or original.get('status')!='ready' or utr.get('status')!='ready':
        raise ValueError('Final annotation requires accepted predictions and UTR validation')
    if (original['gff3_sha256']!=utr['original_model_sha256'] or checksum(prediction/f'{aid}.gff3')!=original['gff3_sha256'] or
        checksum(gff)!=utr['gff3_sha256'] or original['reference_identity']!=utr['reference_identity']):
        raise ValueError('Final annotation source or reference identity differs')
    for field in ('cds_coordinates_and_phases_unchanged','cds_sequences_unchanged','protein_sequences_unchanged'):
        if utr.get(field) is not True:raise ValueError('UTR coding invariance must be established')
    for filename,key in [(f'{aid}_proteins.fa','proteins_sha256'),(f'{aid}_cds.fa','cds_sha256')]:
        if checksum(prediction/filename)!=original[key]:raise ValueError('Original coding sequences changed')
    outdir.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(gff,outdir/f'{aid}.gff3')
    for suffix in ('proteins','cds'):
        shutil.copyfile(prediction/f'{aid}_{suffix}.fa',outdir/f'{aid}_{suffix}.fa')
    status=dict(original,stage='final_annotation',mode=mode,gff3_sha256=utr['gff3_sha256'],
        original_gff3_sha256=original['gff3_sha256'],utr_refinement_skipped=utr['skipped'],
        utr_changed_transcripts=utr['changed_transcripts'],model_set='unfiltered',
        utr_coding_overlap_warning_pairs=utr.get('coding_overlap_warning_pairs',0),
        scope='Validated optional UTR refinement; original coding sequences retained; no model filtering')
    dump(outdir/'annotation_status.json',status)
    (outdir/'README.md').write_text('# Final unfiltered annotation\n\nStatus: **ready**\n\n'
        'These are the accepted models after optional UTR refinement. Original BRAKER outputs remain preserved. '
        'CDS/protein sequences are unchanged. Coding-overlap warnings do not remove models.\n\n'
        f'[Gene models]({aid}.gff3) | [Proteins]({aid}_proteins.fa) | [CDS]({aid}_cds.fa) | '
        '[Validation](annotation_status.json)\n\n'
        f'[UTR review](../../utr/{aid}/README.md) | [Annotation QC](../../../qc/final/README.md) | '
        f'[Native AED](../../aed/{aid}/README.md)\n\n'
        'QC and AED reports are present when those analyses are enabled.\n')
    return status


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','mode','prediction','validation','gff','utr-status','outdir'):p.add_argument('--'+key,required=True)
    a=p.parse_args();finalize(a.assembly_id,a.prediction,a.validation,a.gff,a.utr_status,a.outdir,a.mode)
