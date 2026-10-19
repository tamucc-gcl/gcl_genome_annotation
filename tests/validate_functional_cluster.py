#!/usr/bin/env python3
"""Actual SeqKit/DIAMOND and native merge/GFF decoration; no prediction or downloads."""
import argparse
import csv
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys

from test_functional_annotation import fixture
from validate_utr_cluster import ROOT, execute, dump
from prepare_protein_evidence import checksum
from functional_inputs import TOOLS


def fixtures(root):
    accepted=root/'accepted';dest=accepted/'structural_annotation/identified';rows=[]
    randomizer=random.Random(42)
    seq='M'+''.join(randomizer.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(149))
    for aid,mode in zip('ABCDE',('rna_protein','protein','rna','ab_initio','ab_initio')):
        folder=fixture(dest,aid);protein=folder/f'{aid}_proteins.fa'
        if aid=='E':protein.write_text(''.join(f'>{aid}_g{i}.t1\nMA*AA\n' for i in (1,2,3)))
        else:protein.write_text(f'>{aid}_g1.t1\n{seq}*\n>{aid}_g2.t1\n'+ 'A'*120+f'\n>{aid}_g3.t1\nMA*AA\n')
        assessment=folder/'transcript_assessment.tsv'
        with assessment.open() as handle:
            reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;metrics=list(reader)
        for metric in metrics:
            if mode not in ('rna','rna_protein'):metric['rna_aed']='.'
            if mode not in ('protein','rna_protein'):metric['protein_aed']='.'
        with assessment.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(metrics)
        state=json.loads((folder/'annotation_status.json').read_text());state['proteins_sha256']=checksum(protein)
        state['outputs']['transcript_assessment.tsv']=checksum(assessment);dump(folder/'annotation_status.json',state)
        rows.append(dict(assembly_id=aid,mode=mode,settings={**{f'run_{tool}':tool=='diamond' for tool in TOOLS},'run_omark':False,'eggnog_taxid':'none'}))
    dump(accepted/'preflight/resolved_inputs.json',dict(status='ready',assemblies=rows))
    (accepted/'reference.fa').write_text('>sp|P00001|TEST Test enzyme OS=Synthetic\n'+seq+'\n')


def run_case(root,name,resume=False,fresh=False):
    case=root/('fresh_project' if fresh else 'connected');case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'functional_acceptance.nf'),'-profile','slurm','-w',str(case/'work'),
        '--functional_inputs_root',str(root/'accepted'),'--diamond_database',str(root/'diamond_store'),
        '--diamond_version','synthetic','--n_chunks','2','--outdir',str(case/'results')]
    if resume:command+=['-resume']
    if name=='disabled':command+=['--run_diamond','false']
    if name=='no_hits':command+=['--diamond_evalue','1e-200']
    print(f'Running functional annotation ({name}; actual SeqKit/DIAMOND)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    searches=[t for t in tasks if re.search(r':DIAMOND \(',t['name'])]
    if name=='disabled':assert not searches,f'Disabled DIAMOND unexpectedly ran: {searches}'
    else:
        expected_searches={f'{aid}:eligible.part_{part:03d}' for aid in 'ABCD' for part in (1,2)}
        observed_searches={re.search(r':DIAMOND \((.*)\)$',t['name']).group(1) for t in searches}
        assert len(searches)==8 and observed_searches==expected_searches,(
            f'Expected both protein chunks for each of A/B/C/D (8 DIAMOND tasks); '
            f'observed {len(searches)}: {sorted(observed_searches)}. '
            f'Missing: {sorted(expected_searches-observed_searches)}; trace: {case/"results/pipeline/pipeline_trace.txt"}')
    if name in ('resume','report_change'):
        for t in tasks:
            if name=='report_change' and 'REPORT_FUNCTIONAL_ANNOTATION' in t['name']:continue
            assert t['status']=='CACHED',t
    if name=='fresh_project':
        assert all('DOWNLOAD_DIAMOND_DB' not in t['name'] and 'BUILD_FUNCTIONAL_TEST_DATABASE' not in t['name'] for t in tasks),tasks
    for aid in 'ABCDE':
        folder=case/f'results/functional_annotation/{aid}'
        status=json.loads((folder/'functional_status.json').read_text())
        assert status['status']=='ready' and status['transcripts']==3 and status['coding_structure_unchanged']
        assert status['models_removed']==0
        expected={'disabled':3} if name=='disabled' else {'ineligible':3} if aid=='E' else {'no_hit':2,'ineligible':1} if name=='no_hits' else {'hit':1,'no_hit':1,'ineligible':1}
        assert status['tools']['diamond']==expected,(aid,status['tools']['diamond'],expected)
        with (folder/f'{aid}.functional_annotations.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
        assert len(rows)==3 and all(r['rna_aed']==('0.2' if aid in 'AC' else '.') and
                                  r['protein_aed']==('0.1' if aid in 'AB' else '.') for r in rows)
        assert (folder/'README.md').is_file() and not list(folder.glob('*.html'))
    print('PASS: '+name,flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual SeqKit/DIAMOND; four synthetic accepted evidence modes plus all-ineligible, native shared-store reuse, no-hit/disabled, report/resume cache, parsers for five native formats and unchanged final coding geometry. Full InterProScan/eggNOG/KofamScan/DeepSig database pilot pending.')
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
        for name in ('main.nf','functional_annotation.nf','nextflow.config'):shutil.copy2(ROOT/name,pipeline/name)
        shutil.copy2(ROOT/'tests/functional_fixture.nf',pipeline/'functional_acceptance.nf')
        with (pipeline/'nextflow.config').open('a') as handle:
            handle.write("\nprocess {\n withLabel: diamond { cpus=2; memory='4 GB'; time='1h' }\n withLabel: chunk_proteins { cpus=1; memory='2 GB'; time='1h' }\n}\n")
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        text=probe.stdout+probe.stderr;(root/'main_syntax.console.log').write_text(text)
        assert probe.returncode and 'Provide --sample_sheet' in text,text[-8000:]
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        helper=pipeline/'py_scripts/functional_report.py';helper.write_text(helper.read_text()+'\n# report-only cache fixture\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'fresh_project',fresh=True))
        summary['checks'].append(run_case(root,'no_hits',True))
        summary['checks'].append(run_case(root,'disabled',True))
        summary['status']='passed'
    except Exception as exc:
        summary['error']=str(exc);print('FAIL: '+str(exc),flush=True)
    dump(root/'acceptance.json',summary);print('Acceptance summary: '+str(root/'acceptance.json'))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
