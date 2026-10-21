#!/usr/bin/env python3
"""Actual Markdown/Matplotlib report using accepted synthetic releases; no searches."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

from validate_utr_cluster import ROOT, dump, execute
from validate_release_cluster import fixtures as release_fixtures
from test_unified_report import report_fixture
from reference_repeats import sha


def run_case(root,name,resume):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'report_acceptance.nf'),'-profile','slurm','-w',str(case/'work'),
        '--report_test_request',str(root/'report_request.json'),'--outdir',str(case/'results')]
    if resume:command+=['-resume']
    print(f'Running unified Markdown report ({name}; no searches)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    assert len(tasks)==1 and tasks[0]['name']=='ANNOTATION_REPORT',tasks
    if name=='resume':assert tasks[0]['status']=='CACHED',tasks
    if name=='report_change':assert tasks[0]['status']=='COMPLETED',tasks
    folder=case/'results';state=json.loads((folder/'reports/annotation_report_status.json').read_text())
    text=(folder/'annotation_report.md').read_text()
    assert state['status']=='ready' and state['assemblies']==['A','B'],state
    assert '## Assembly 1: A' in text and '## Assembly 2: B' in text
    assert 'Tuning Recommendations' not in text and '<html' not in text
    for file,digest in state['outputs'].items():
        path=folder/file if file=='annotation_report.md' else folder/'reports'/file
        assert sha(path)==digest,file
    # All generated plot links must resolve; native structural images remain
    # external stage links, as in the complete output hierarchy.
    for link in re.findall(r'\]\((reports/annotation_report_plots/[^)]+)\)',text):assert (folder/link).is_file(),link
    assert not list((folder/'reports').rglob('*.html'))
    print('PASS: '+name,flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty or use --resume')
    summary=dict(status='failed',checks=[],scope='Actual Nextflow/Matplotlib unified Markdown report; two accepted synthetic releases with mixed RNA/protein and genome-only modes, missing evidence/disabled tracks, current manifest counts, source checksums, PNG/SVG, unchanged resume and report-only rerun. No scientific searches; full-data report pilot pending.')
    try:
        execute([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],ROOT,root/'unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume:
            release_fixtures(root)
            sources=json.loads((root/'requests.json').read_text())
            request=report_fixture(root,mode=dict(A='rna_protein',B='ab_initio'),sources=sources)
            dump(root/'report_request.json',request)
        pipeline=root/'pipeline';pipeline.mkdir(exist_ok=True)
        for directory in ('modules','workflows','py_scripts'):
            for source in (ROOT/directory).rglob('*'):
                if not source.is_file() or '__pycache__' in source.parts:continue
                dest=pipeline/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists() or source.read_bytes()!=dest.read_bytes():shutil.copy2(source,dest)
        shutil.copy2(ROOT/'main.nf',pipeline/'main.nf');shutil.copy2(ROOT/'nextflow.config',pipeline/'nextflow.config')
        shutil.copy2(ROOT/'tests/annotation_report_fixture.nf',pipeline/'report_acceptance.nf')
        with (pipeline/'nextflow.config').open('a') as h:h.write('\nparams.report_test_request=null\n')
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        text=probe.stdout+probe.stderr;(root/'main_syntax.console.log').write_text(text)
        assert probe.returncode and 'Provide --sample_sheet' in text,text[-5000:]
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        helper=pipeline/'py_scripts/unified_report.py';helper.write_text(helper.read_text().replace('Results are reported per assembly.','Results are reported per assembly. '))
        summary['checks'].append(run_case(root,'report_change',True))
        summary['status']='passed'
    except Exception as exc:
        summary.update(error=str(exc),traceback=traceback.format_exc());print('FAIL: '+str(exc),flush=True)
    dump(root/'acceptance.json',summary);print('Acceptance summary: '+str(root/'acceptance.json'))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
