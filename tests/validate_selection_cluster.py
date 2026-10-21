#!/usr/bin/env python3
"""Nextflow model selection on tiny accepted fixtures; no evidence recalculation."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from test_model_selection import fixture, change_assessment
from validate_utr_cluster import execute, dump, ROOT
sys.path.insert(0,str(ROOT/'py_scripts'))
from prepare_protein_evidence import checksum

MODES=dict(A='rna_protein',B='protein',C='rna',D='ab_initio')


def fixtures(root):
    accepted=root/'accepted'; predictions=[]
    for aid,mode in MODES.items():
        source=root/f'source_{aid}';source.mkdir()
        g,p,c,identity,validation,assessment=fixture(source,single=True,shared=True)
        state=json.loads(validation.read_text());state.update(assembly_id=aid,mode=mode)
        validation.write_text(json.dumps(state))
        folder=accepted/f'structural_annotation/final/{aid}';folder.mkdir(parents=True)
        for src,name in [(g,f'{aid}.gff3'),(p,f'{aid}_proteins.fa'),(c,f'{aid}_cds.fa'),(validation,'annotation_status.json')]:
            shutil.copy2(src,folder/name)
        target=accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json';target.parent.mkdir(parents=True);shutil.copy2(identity,target)
        rna=mode in ('rna','rna_protein');protein=mode in ('protein','rna_protein')
        changes={tid:dict(rna_aed=str(value) if rna else '.',protein_aed=str(value) if protein else '.',
            te_gene_flag=str(tid!='t3'),junction_coding_junctions='0' if tid=='t3' else '1',
            junction_coding_pass_fraction='.' if tid=='t3' or not rna else ('1.0' if tid=='t1' else '0.5')) for tid,value in [('t1',.1),('t2',.8),('t3',.1)]}
        change_assessment(assessment,changes,aed_enabled=True,junction_qc_enabled=True,te_enabled=True,assembly_id=aid)
        target=accepted/f'structural_annotation/assessment/{aid}';shutil.copytree(assessment,target)
        predictions.append(dict(assembly_id=aid,mode=mode,status='ready'))
    dump(accepted/'structural_annotation/prediction_status.json',dict(status='ready',assemblies=predictions))


def run_case(root,name,resume=False):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'validation/entrypoints/model_selection.nf'),'-profile','slurm','-w',str(case/'work'),
        '--selection_inputs_root',str(root/'accepted'),'--outdir',str(case/'results')]
    if resume:command.append('-resume')
    if name=='threshold_change':command+=['--filter_min_cds_length','50']
    if name in ('aed_rna','aed_protein'):
        command+=['--filter_structural_action','off','--filter_max_aed','5e-1','--filter_aed_source',name[4:]]
    if name=='junction':command+=['--filter_structural_action','off','--filter_min_junction_fraction','0.9']
    if name=='te_remove':command+=['--filter_structural_action','off','--filter_te_action','remove']
    if name=='flag_only':command+=['--filter_structural_action','flag']
    if name=='disabled':command+=['--run_model_selection','false']
    print(f'Running model selection ({name}; no scientific searches)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    trace=case/'results/pipeline/pipeline_trace.txt'
    with trace.open() as handle:tasks=list(csv.DictReader(handle,delimiter='\t'))
    jobs=[t for t in tasks if 'FILTER_GENE_MODELS (' in t['name']]
    if name=='disabled':
        assert not jobs,jobs
        print('PASS: disabled',flush=True);return dict(case=name,passed=True)
    assert len(jobs)==4 and len(tasks)==4,tasks
    if name=='resume':assert all(t['status']=='CACHED' for t in jobs),jobs
    for aid,mode in MODES.items():
        folder=case/f'results/structural_annotation/selection/{aid}'
        state=json.loads((folder/'selection_status.json').read_text())
        assert state['status']=='ready' and state['input_transcripts']==3 and state['input_genes']==2
        expected={'t1','t2'}
        if name in ('threshold_change','flag_only'):expected={'t1','t2','t3'}
        if name in ('aed_rna','aed_protein','junction'):
            available=(mode in ('rna','rna_protein')) if name!='aed_protein' else (mode in ('protein','rna_protein'))
            expected={'t1','t3'} if available else {'t1','t2','t3'}
        if name=='te_remove':expected={'t3'}
        with (folder/'retained/transcript_assessment.tsv').open() as handle:retained=list(csv.DictReader(handle,delimiter='\t'))
        assert {r['transcript_id'] for r in retained}==expected,(aid,name,retained)
        assert state['retained_transcripts']==len(expected)
        assert state['excluded_transcripts']==3-len(expected)
        for filename,sha in state['outputs'].items():assert checksum(folder/filename)==sha,filename
        for label in ('retained','excluded'):
            text=(folder/f'{label}/{aid}.gff3').read_text()
            if name in ('aed_rna','aed_protein','junction') and len(expected)==2:
                assert 'Parent=t1,t2' not in text
    print(f'PASS: {name}',flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual Nextflow selection; synthetic accepted four-mode assessment handoffs, original defaults, source-specific AED, coding junction/TE opt-in, shared parents, synchronized outputs, resume/threshold cache; no mapping, repeats, BRAKER, AED or junction recalculation')
    try:
        version=subprocess.run(['nextflow','-version'],capture_output=True,text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b',version.stdout+version.stderr):raise RuntimeError('Load Nextflow 23.10.1')
        execute([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],ROOT,root/'unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume:fixtures(root)
        pipeline=root/'pipeline';pipeline.mkdir(exist_ok=True)
        for directory in ('modules','workflows','py_scripts','environments'):
            for source in (ROOT/directory).rglob('*'):
                if not source.is_file() or '__pycache__' in source.parts:continue
                dest=pipeline/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists() or source.read_bytes()!=dest.read_bytes():shutil.copy2(source,dest)
        for name in ('main.nf','validation/entrypoints/model_selection.nf','nextflow.config'):
            (pipeline/name).parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(ROOT/name,pipeline/name)
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        (root/'main_syntax.console.log').write_text(probe.stdout+probe.stderr)
        if probe.returncode==0 or 'Provide --sample_sheet' not in probe.stdout+probe.stderr:raise RuntimeError('main syntax/input guard failed\n'+(probe.stdout+probe.stderr)[-6000:])
        summary['checks'].append(dict(case='main_syntax',passed=True))
        for name in ['recover' if a.resume else 'fresh','resume','threshold_change','aed_rna','aed_protein','junction','te_remove','flag_only','disabled']:
            summary['checks'].append(run_case(root,name,resume=(a.resume or name!='fresh')))
        summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root / "acceptance.json"}',flush=True)
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
