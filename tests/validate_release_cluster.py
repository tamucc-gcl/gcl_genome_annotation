#!/usr/bin/env python3
"""Actual Nextflow packaging of accepted synthetic inputs; no scientific searches."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

from validate_utr_cluster import ROOT, dump, execute
from test_annotation_release import release_fixture
from reference_repeats import sha


def fixtures(root):
    rows=[]
    for aid in ('A','B'):
        folder=root/'fixtures'/aid; folder.mkdir(parents=True)
        request=release_fixture(folder)
        if aid=='B':
            # A second independent, accepted model set, without optional tracks.
            identified=Path(request['identified'])
            for p in list(identified.iterdir()):
                if p.is_file() and p.name!='annotation_status.json':
                    text=p.read_text().replace('A_g','B_g')
                    p.write_text(text)
                    if p.name.startswith('A.') or p.name.startswith('A_'):p.rename(p.with_name('B'+p.name[1:]))
            state=json.loads((identified/'annotation_status.json').read_text())
            state.update(assembly_id='B',prefix='B',outputs={p.name:sha(p) for p in identified.iterdir() if p.is_file() and p.name not in ('annotation_status.json','README.md')})
            for suffix,key in [('gff3','gff3_sha256'),('proteins.fa','proteins_sha256'),('cds.fa','cds_sha256')]:
                state[key]=sha(identified/('B.gff3' if suffix=='gff3' else 'B_'+suffix))
            dump(identified/'annotation_status.json',state)
            for key in ('track_status','repeat_status'):
                p=Path(request[key]); state=json.loads(p.read_text()); state['assembly_id']='B'
                if key=='track_status':state['tracks']={t:dict(status='disabled') for t in state['tracks']}
                dump(p,state)
            request.update(aid='B',tracks=[])
            for key in ('gff','table','genes'):
                p=Path(request['functional'][key]);p.write_text(p.read_text().replace('A_g','B_g'))
            f=request['functional'];p=Path(f['status']);state=json.loads(p.read_text())
            state.update(assembly_id='B',original_gff3_sha256=sha(identified/'B.gff3'),annotated_gff3_sha256=sha(f['gff']),table_sha256=sha(f['table']))
            dump(p,state)
        rows.append(request)
    dump(root/'requests.json',rows)


def run_case(root,name,resume):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'release_acceptance.nf'),'-profile','slurm','-w',str(case/'work'),
        '--release_test_requests',str(root/'requests.json'),'--outdir',str(case/'results')]
    if resume:command+=['-resume']
    if name=='functional_disabled':command+=['--run_functional_annotation','false']
    print(f'Running final packaging ({name}; no searches)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    packages=[t for t in tasks if t['name'].split(' (')[0].endswith('PACKAGE_ANNOTATION_RELEASE')]
    assert len(packages)==2,packages
    if name=='resume':assert all(t['status']=='CACHED' for t in tasks),tasks
    if name=='index_change':assert all(t['status']=='CACHED' for t in packages),packages
    release=case/'results/release'
    for aid in ('A','B'):
        folder=release/aid;state=json.loads((folder/'release_status.json').read_text())
        assert state['status']=='ready' and state['assembly_id']==aid,state
        assert state['functional_enabled']==(name!='functional_disabled'),state
        for file,metadata in state['outputs'].items():assert sha(folder/file)==metadata['sha256'],file
        assert state['disabled_tracks']==(['centromere'] if aid=='A' else ['trnascan','infernal','telomere','centromere']),state
        assert not list(folder.glob('*.html'))
    index=json.loads((release/'release_index.json').read_text())
    assert [r['assembly_id'] for r in index['assemblies']]==['A','B'],index
    print('PASS: '+name,flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty or use --resume')
    summary=dict(status='failed',checks=[],scope='Actual Nextflow release packaging; two accepted synthetic assemblies, selected/disabled/empty tracks, functional on/off, complete coding/track ID maps, checksums, parent/reference checks, resume and index-only cache. No mapping, repeats, prediction or functional searches; full-genome release pilot pending.')
    try:
        execute([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],ROOT,root/'unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume:fixtures(root)
        pipeline=root/'pipeline';pipeline.mkdir(exist_ok=True)
        for directory in ('modules','workflows','py_scripts'):
            for source in (ROOT/directory).rglob('*'):
                if not source.is_file() or '__pycache__' in source.parts:continue
                dest=pipeline/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists() or source.read_bytes()!=dest.read_bytes():shutil.copy2(source,dest)
        shutil.copy2(ROOT/'main.nf',pipeline/'main.nf');shutil.copy2(ROOT/'nextflow.config',pipeline/'nextflow.config')
        shutil.copy2(ROOT/'tests/release_fixture.nf',pipeline/'release_acceptance.nf')
        with (pipeline/'nextflow.config').open('a') as h:h.write('\nparams.release_test_requests=null\n')
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        text=probe.stdout+probe.stderr;(root/'main_syntax.console.log').write_text(text)
        assert probe.returncode and 'Provide --sample_sheet' in text,text[-5000:]
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        module=pipeline/'modules/annotationRelease.nf';module.write_text(module.read_text().replace('# Final annotation releases','# Final annotation releases (index cache test)'))
        summary['checks'].append(run_case(root,'index_change',True))
        summary['checks'].append(run_case(root,'functional_disabled',True))
        summary['checks'].append(run_case(root,'reenabled',True))
        summary['status']='passed'
    except Exception as exc:
        summary.update(error=str(exc),traceback=traceback.format_exc());print('FAIL: '+str(exc),flush=True)
    dump(root/'acceptance.json',summary);print('Acceptance summary: '+str(root/'acceptance.json'))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
