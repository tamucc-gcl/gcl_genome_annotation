#!/usr/bin/env python3
"""Native tRNAscan/Infernal/tidk/centroAnno acceptance without prediction or repeats."""
import argparse
import csv
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys

from validate_utr_cluster import ROOT, dump, execute
from reference_repeats import reference


def task_process(task):
    """Match the process itself, excluding workflow scopes and the task tag."""
    return task['name'].split(' (', 1)[0].rsplit(':', 1)[-1]


def check_telomere_ends(folder, aid, identity, motif):
    """Both inserted terminal arrays must be covered; a continuous call can span both."""
    lengths={c['id']:c['length'] for c in json.loads(Path(identity).read_text())['contigs']}
    rows=[line.split('\t') for line in (folder/f'{aid}.telomere.gff3').read_text().splitlines()
          if line.strip() and not line.startswith('#')]
    terminal=[r for r in rows if r[0]=='chr1' and r[2]=='telomeric_repeat']
    assert any(int(r[3])==1 for r in terminal),f'{aid}: missing left terminal array'
    assert any(int(r[4])==lengths['chr1'] for r in terminal),f'{aid}: missing right terminal array'
    for row in rows:
        attrs=dict(token.split('=',1) for token in row[8].split(';') if '=' in token)
        assert attrs['motif']==motif,f'{aid}: wrong telomere motif: {attrs}'


def check_centromere_mode(log):
    # Native package output: the default satellite-only mode never performs this scan.
    text=Path(log).read_text()
    assert 'Detecting tandem repeat regions in the ' in text,(
        f'{log}: native assembly-region scan not found; verify packaged centroAnno uses -x anno-asm')


def fixtures(root):
    folder=root/'fixtures';folder.mkdir(exist_ok=True)
    rng=random.Random(42)
    # No ambiguous homopolymer input: provide tandem satellite arrays and a known model sequence.
    unit=''.join(rng.choice('ACGT') for _ in range(171))
    rna='GCGCGC'+'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA'+'GCGCGC'
    rows=[]
    for aid,motif in (('A','TTAGGG'),('B','TTAGG'),('C','TTAGGG'),('D','TTAGGG')):
        source=folder/f'{aid}.source.fa'
        source.write_text('>chr1\n'+motif*60+unit*150+rna.replace('U','T')+motif*60+'\n>short\n'+'AC'*100+'\n')
        if aid=='D':source.write_text('>chr1\n'+'AC'*100+'\n')
        bundle=folder/aid
        state=reference(source,aid,0,bundle)
        assert state['status']=='ready'
        shutil.copy2(bundle/f'{aid}.filtered.fasta',folder/f'{aid}.fa')
        shutil.copy2(bundle/f'{aid}.filtered.fasta',folder/f'{aid}.masked.fa')
        shutil.copy2(bundle/f'{aid}.contigs.json',folder/f'{aid}.json')
        rows.append(dict(assembly_id=aid,settings={**{'run_'+t:aid in 'AB' or (aid=='D' and t=='telomere') for t in ('trnascan','infernal','telomere','centromere')},'telomere_motif':motif}))
    dump(folder/'assemblies.json',rows)
    (folder/'model.sto').write_text('# STOCKHOLM 1.0\n#=GF ID Test_RNA\n#=GF AC RF99999\n#=GF GA 5.0\ns1 '+rna+'\n#=GC SS_cons '+'('*6+'.'*(len(rna)-12)+')'*6+'\n//\n')


def run_case(root,name,resume=False):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'genome_tracks_acceptance.nf'),'-profile','slurm','-w',str(case/'work'),
        '--genome_tracks_test_root',str(root/'fixtures'),'--rfam_database',str(root/'rfam_store'),
        '--rfam_version','synthetic','--n_infernal_chunks','3','--telomere_window','1000',
        '--telomere_end_window','1000','--outdir',str(case/'results')]
    if resume:command+=['-resume']
    if name=='disabled':command+=['--genome_tracks_test_disabled','true']
    print(f'Running genome tracks ({name}; actual tRNAscan/Infernal/tidk/centroAnno)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    scientific=('tRNASCAN','INFERNAL','TELOMERE_ANNOTATION','CENTROMERE_ANNOTATION')
    if name=='disabled':assert not any(task_process(r) in scientific for r in tasks)
    else:
        for tool in scientific:
            matches=[r for r in tasks if task_process(r)==tool]
            assert len(matches)==(4 if tool=='INFERNAL' else 3 if tool=='TELOMERE_ANNOTATION' else 2),(tool,matches)
        if name in ('resume','report_change'):
            for task in tasks:
                if name=='report_change' and 'REPORT_GENOME_TRACKS' in task['name']:continue
                assert task['status']=='CACHED',task
    for aid in 'ABCD':
        folder=case/f'results/genome_tracks/{aid}'
        state=json.loads((folder/'genome_tracks_status.json').read_text())
        assert state['status']=='ready'
        disabled=name=='disabled' or aid=='C'
        assert all(s['status']==('disabled' if disabled or (aid=='D' and t!='telomere') else 'ready') for t,s in state['tracks'].items())
        if aid=='D' and not disabled:assert state['tracks']['telomere']['features']==0,state
        if not disabled and aid!='D':
            check_centromere_mode(case/f'results/centromere/{aid}.centroanno.log')
            assert state['tracks']['telomere']['features']>=1,state
            check_telomere_ends(folder,aid,root/'fixtures'/f'{aid}.json',state['settings']['telomere_motif'])
            assert state['tracks']['infernal']['features']>=1,state
            for track in state['tracks']:
                assert (folder/f'{aid}.{track}.gff3').is_file()
        assert not list(folder.glob('*.html'))
    print('PASS: '+name,flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual tRNAscan-SE, synthetic calibrated Infernal CM, tidk and centroAnno 1.0.2; two assemblies/motifs, fewer contigs than requested chunks, disabled/no-hit tracks, coordinates/IDs, resume/report-only cache; no BRAKER, mapping, repeats or functional searches; full Rfam/genome pilot pending.')
    try:
        execute([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],ROOT,root/'unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume:fixtures(root)
        pipeline=root/'pipeline';pipeline.mkdir(exist_ok=True)
        for directory in ('modules','workflows','py_scripts','perl_scripts'):
            for source in (ROOT/directory).rglob('*'):
                if not source.is_file() or '__pycache__' in source.parts:continue
                dest=pipeline/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists() or source.read_bytes()!=dest.read_bytes():shutil.copy2(source,dest)
        shutil.copy2(ROOT/'main.nf',pipeline/'main.nf');shutil.copy2(ROOT/'nextflow.config',pipeline/'nextflow.config')
        shutil.copy2(ROOT/'tests/genome_tracks_fixture.nf',pipeline/'genome_tracks_acceptance.nf')
        with (pipeline/'nextflow.config').open('a') as h:
            h.write("\nparams.genome_tracks_test_disabled=false\nprocess {\n withLabel: download_rfam_db { cpus=2; memory='4 GB'; time='2h' }\n withLabel: infernal { cpus=2; memory='4 GB'; time='2h' }\n withLabel: trna_scan { cpus=2; memory='4 GB'; time='2h' }\n withLabel: centromere { cpus=2; memory='8 GB'; time='2h' }\n withName: BUILD_TRACK_TEST_RFAM { scratch=false }\n}\n")
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        text=probe.stdout+probe.stderr;(root/'main_syntax.console.log').write_text(text)
        assert probe.returncode and 'Provide --sample_sheet' in text,text[-8000:]
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        # Only the report helper changes: expensive native tools must remain cached.
        helper=pipeline/'py_scripts/report_genome_tracks.py';helper.write_text(helper.read_text()+'\n# report-only cache fixture\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'disabled',True))
        summary['status']='passed'
    except Exception as exc:
        summary['error']=str(exc);print('FAIL: '+str(exc),flush=True)
    dump(root/'acceptance.json',summary);print('Acceptance summary: '+str(root/'acceptance.json'))
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
