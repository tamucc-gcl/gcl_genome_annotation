#!/usr/bin/env python3
"""Actual UTR/final-QC/native-AED integration on tiny accepted model fixtures."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from validate_utr_cluster import fixtures as utr_fixtures, execute, dump, ROOT
sys.path.insert(0,str(ROOT/'py_scripts'))
from prepare_protein_evidence import checksum


def fixtures(root):
    utr_fixtures(root)
    accepted=root/'accepted';rows=[]
    for aid,strand in [('A','+'),('C','-')]:
        sam=root/f'sams/{aid}.sam'
        sam.write_text(''.join(line if line.startswith('@') else line.rstrip('\n')+f'\tNH:i:1\tXS:A:{strand}\n'
            for line in sam.read_text().splitlines(keepends=True)))
    # RNA-only BRAKER does not retain ETP StringTie output: exercise assembly fallback.
    (accepted/'structural_annotation/braker/C/braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff').unlink()
    for aid in ('A','B','C','D'):
        if aid=='D':
            for source,target in [
                (accepted/'structural_annotation/braker/B/B.gff3',accepted/'structural_annotation/braker/D/D.gff3'),
                (accepted/'repeats/B.masked.fasta',accepted/'repeats/D.masked.fasta'),
                (accepted/'filtered_assemblies/B/B.contigs.json',accepted/'filtered_assemblies/D/D.contigs.json')]:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
            row=dict(assembly_id=aid,status='ready',mode='ab_initio',readsets=[])
        else:row=json.loads((accepted/f'structural_annotation/braker_status/{aid}.prediction.json').read_text())
        directory=accepted/f'structural_annotation/braker/{aid}'
        amino='F' if aid=='C' else 'K';nucleotide='T' if aid=='C' else 'A'
        proteins=directory/f'{aid}_proteins.fa';proteins.write_text('>t1\n'+amino*100+'\n')
        cds=directory/f'{aid}_cds.fa';cds.write_text('>t1\n'+nucleotide*300+'\n')
        reference=json.loads((accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json').read_text())
        row.update(genes=1,coding_transcripts=1,reference_identity=reference['reference_identity'],
            gff3_sha256=checksum(directory/f'{aid}.gff3'),proteins_sha256=checksum(proteins),cds_sha256=checksum(cds))
        dump(accepted/f'structural_annotation/braker_status/{aid}.prediction.json',row)
        rows.append(row)
        if row['mode'] in ('rna_protein','protein'):
            source=accepted/f'protein_evidence/{aid}/proteins.fa';source.parent.mkdir(parents=True)
            source.write_text('>evidence1\n'+'K'*100+'\n')
            dump(source.parent/'evidence.json',dict(status='ready',proteins_sha256=checksum(source)))
    dump(accepted/'fixtures.json',rows)


def run_case(root,name,resume=False,enabled=True):
    pipeline=root/'pipeline';case=root/'connected';case.mkdir(exist_ok=True)
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'postprocessing_fixture.nf'),'-profile','slurm','-w',str(case/'work'),
        '--utr_inputs_root',str(root/'accepted'),'--outdir',str(case/'results'),
        '--utr_refine',str(enabled).lower(),'--utr_ranking_cpus','2','--utr_ranking_memory','4 GB',
        '--aed_alignment_cpus','2','--aed_alignment_memory','4 GB',
        '--ingenannot_aed_cpus','2','--ingenannot_aed_memory','4 GB']
    if resume:command.append('-resume')
    print(f'Running connected UTR/final QC/native AED ({name})...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    for aid,mode in [('A','rna_protein'),('B','protein'),('C','rna'),('D','ab_initio')]:
        final=case/f'results/structural_annotation/final/{aid}'
        state=json.loads((final/'annotation_status.json').read_text())
        assert state['status']=='ready' and state['mode']==mode
        assert state['utr_refinement_skipped']==(not enabled or aid in ('B','D'))
        assert state['utr_changed_transcripts']==(1 if enabled and aid in ('A','C') else 0)
        assert state['gff3_sha256']==checksum(final/f'{aid}.gff3')
        assert checksum(final/f'{aid}_proteins.fa')==checksum(root/f'accepted/structural_annotation/braker/{aid}/{aid}_proteins.fa')
        scoring=case/f'results/structural_annotation/aed/{aid}'
        assessed=json.loads((scoring/'aed_status.json').read_text())
        assert assessed['status']=='ready' and assessed['model_sha256']==state['gff3_sha256']
        assert (scoring/'native/native.aed.gff3').is_file()
        with (scoring/'aed_scores.tsv').open() as h:row=next(csv.DictReader(h,delimiter='\t'))
        assert (row['rna_aed']!='.')==(mode in ('rna_protein','rna'))
        assert (row['protein_aed']!='.')==(mode in ('rna_protein','protein'))
        assert not any('combined' in key or 'legacy' in key for key in row)
        selection=json.loads((case/f'results/qc/final/{aid}/selection/selection.json').read_text())
        assert selection['gff3_sha256']==state['gff3_sha256'] and selection['model_source']=='final_annotation'
    qc=json.loads((case/'results/qc/final/baseline_status.json').read_text())
    assert qc['status']=='ready' and len(qc['assemblies'])==4
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    heavy=[t for t in tasks if re.search(r':(RUN_UTR_REFINEMENT|RUN_INGENANNOT_AED|ALIGN_AED_PROTEINS) \(',t['name'])]
    assert len(heavy)==10,f'Expected 4 UTR, 4 AED, 2 protein tasks: {heavy}'
    if name in ('resume','report_change'):
        assert all(t['status']=='CACHED' for t in heavy),f'Unexpected cache miss: {heavy}'
    if name=='disabled':
        assert all(t['status']=='CACHED' for t in heavy if ':ALIGN_AED_PROTEINS ' in t['name'])
    assert not any(re.search(r'(REPEAT_MODEL|ALIGN_RNA_READSET|:BRAKER \()',t['name']) for t in tasks)
    print(f'PASS: connected annotation {name}',flush=True)
    return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual tiny UTR/gffread/AGAT/miniprot/native AED; all four evidence modes, main syntax, resume and reporting cache; no BRAKER training or repeats')
    try:
        version=subprocess.run(['nextflow','-version'],capture_output=True,text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b',version.stdout+version.stderr):raise RuntimeError('Load Nextflow 23.10.1')
        with (root/'unit_tests.log').open('w') as log:
            result=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError('Python tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume:fixtures(root)
        pipeline=root/'pipeline';pipeline.mkdir(exist_ok=True)
        for directory in ('modules','workflows','py_scripts','environments'):
            for source in (ROOT/directory).rglob('*'):
                if not source.is_file() or '__pycache__' in source.parts:continue
                dest=pipeline/source.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists() or source.read_bytes()!=dest.read_bytes():shutil.copy2(source,dest)
        for name in ('main.nf','nextflow.config'):shutil.copy2(ROOT/name,pipeline/name)
        shutil.copy2(ROOT/'tests/postprocessing_fixture.nf',pipeline/'postprocessing_fixture.nf')
        # No sample sheet: main must compile and reach its explicit input guard, without submitting work.
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        (root/'main_syntax.console.log').write_text(probe.stdout+probe.stderr)
        if probe.returncode==0 or 'Provide --sample_sheet' not in probe.stdout+probe.stderr:raise RuntimeError('main.nf syntax/input guard failed; inspect main_syntax.console.log')
        summary['checks'].append(dict(case='main_syntax',passed=True))
        execute(['nextflow','-C',str(ROOT/'nextflow.config'),'run',str(ROOT/'tests/utr_fixture_bams.nf'),'-profile','slurm',
            '--fixture_sams',str(root/'sams'),'--outdir',str(root/'accepted'),'-w',str(root/'bam_work'),
            '--rna_alignment_cpus','2','--rna_alignment_memory','4 GB','-resume'],root,root/'bam_fixture.console.log')
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        with (pipeline/'py_scripts/final_annotation.py').open('a') as h:h.write('\n# Reporting cache regression.\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'disabled',True,False))
        summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root/"acceptance.json"}')
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
