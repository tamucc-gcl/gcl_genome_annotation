#!/usr/bin/env python3
"""Compare completed native AED output with accepted models without rescoring."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from evidence_aed import features, models, merged
from prepare_protein_evidence import checksum
from utr_refinement import subtract


def difference(before,after):
    if sorted(before)==sorted(after):return 'identical'
    if set(before)==set(after):return 'duplicate_records_only'
    if merged(before)==merged(after):return 'same_bases_different_segmentation'
    return 'changed_bases'


def inspect(original,native,identity,validation,output,focus=()):
    reference=json.loads(Path(identity).read_text());accepted=json.loads(Path(validation).read_text())
    if accepted['status']!='ready' or checksum(original)!=accepted['gff3_sha256'] or reference['reference_identity']!=accepted['reference_identity']:
        raise ValueError('Diagnostic inputs differ from accepted final annotation')
    lengths={c['id']:c['length'] for c in reference['contigs']}
    print('Reading accepted and completed native models...',flush=True)
    before,after=models(original,lengths),models(native,lengths)
    def phases(path):
        result={}
        for row in features(path,lengths):
            if row['kind']=='CDS':
                for tid in row['attrs']['Parent'].split(','):
                    result.setdefault(tid,[]).append((row['start'],row['end'],row['phase']))
        return result
    old_phases,new_phases=phases(original),phases(native)
    categories=Counter();details=[];ids=[]
    for tid in sorted(set(before)&set(after)):
        a,b=before[tid],after[tid]
        changes=[]
        for key in ('gene_id','seq','strand'):
            if a[key]!=b[key]:changes.append(key+'_changed')
        for key in ('exons','cds'):
            category=difference(a[key],b[key])
            if category!='identical':changes.append(key+':'+category)
        if sorted(old_phases.get(tid,[]))!=sorted(new_phases.get(tid,[])):
            category='duplicate_records_only' if set(old_phases.get(tid,[]))==set(new_phases.get(tid,[])) else 'changed'
            changes.append('cds_phases:'+category)
        bounds=(a['start'],a['end'])!=(b['start'],b['end'])
        if bounds and (b['start'],b['end'])!=(min(s for s,e in a['exons']),max(e for s,e in a['exons'])):
            changes.append('unexpected_transcript_bounds')
        if changes:
            categories.update(changes);ids.append(tid)
        if (changes and len(details)<10) or tid in focus:
            details.append(dict(transcript_id=tid,changes=changes,
                original={k:a[k] for k in ('gene_id','seq','strand','start','end','exons','cds','attrs')},
                native={k:b[k] for k in ('gene_id','seq','strand','start','end','exons','cds','attrs')},
                original_cds_phases=old_phases.get(tid,[]),native_cds_phases=new_phases.get(tid,[]),
                added_exonic_bases=subtract(b['exons'],a['exons']),removed_exonic_bases=subtract(a['exons'],b['exons'])))
    result=dict(stage='native_aed_structure_diagnostic',coordinates='0-based half-open',
        scope='Read-only completed-output comparison; no scoring or annotation changes',
        original_model_sha256=checksum(original),native_output_sha256=checksum(native),
        original_transcripts=len(before),native_transcripts=len(after),
        missing_native_transcripts=sorted(set(before)-set(after)),extra_native_transcripts=sorted(set(after)-set(before)),
        mismatched_transcripts=len(ids),categories=dict(categories),mismatched_transcript_ids=ids,examples=details)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('examples','mismatched_transcript_ids')},indent=2))
    print(f'Diagnostic report: {output}',flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs-root',required=True);p.add_argument('--assembly-id',required=True)
    p.add_argument('--output',required=True);p.add_argument('--transcript-id',action='append',default=[])
    a=p.parse_args();root=Path(a.inputs_root);aid=a.assembly_id
    inspect(root/f'structural_annotation/final/{aid}/{aid}.gff3',
        root/f'structural_annotation/aed/{aid}/native/native.aed.gff3',
        root/f'filtered_assemblies/{aid}/{aid}.contigs.json',
        root/f'structural_annotation/final/{aid}/annotation_status.json',a.output,a.transcript_id)
