#!/usr/bin/env python3
"""Actual Portcullis fixture: strand/coordinates, no RNA, cache and model-only change."""
import argparse
import csv
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from prepare_protein_evidence import checksum
from reference_repeats import reference


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')


def execute(command,cwd,log):
    with log.open('w') as h:result=subprocess.run(command,cwd=cwd,stdout=h,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'{log}: exit {result.returncode}\n'+log.read_text(errors='replace')[-6000:])


def fixtures(root):
    rng=random.Random(64);dna=list(''.join(rng.choice('ACGT') for _ in range(2000)))
    for start,end,strand in [(200,300,'+'),(400,500,'+'),(1100,1200,'-')]:
        dna[start:start+2]=list('GT' if strand=='+' else 'CT');dna[end-2:end]=list('AG' if strand=='+' else 'AC')
    dna=''.join(dna);accepted=root/'accepted';rows=[]
    model='##gff-version 3\n'
    for tid,strand,exons,cds in [('t1','+',[(100,200),(300,400),(500,650)],[(150,200),(300,350)]),
                               ('t2','-',[(1000,1100),(1200,1400)],[(1000,1100),(1200,1400)]),
                               ('t3','+',[(1600,1800)],[(1600,1800)])]:
        for kind,s,e,phase,attrs in [('gene',exons[0][0],exons[-1][1],'.','ID=g'+tid),
                                   ('mRNA',exons[0][0],exons[-1][1],'.',f'ID={tid};Parent=g{tid}')]+[
            (kind,s,e,'0' if kind=='CDS' else '.',f'ID={tid}.{kind}{i};Parent={tid}')
            for kind,parts in [('exon',exons),('CDS',cds)] for i,(s,e) in enumerate(parts)]:
            model+=f'chr1\tfixture\t{kind}\t{s+1}\t{e}\t.\t{strand}\t{phase}\t{attrs}\n'
    for aid,mode in [('A','rna_protein'),('B','ab_initio')]:
        genome=accepted/f'repeats/{aid}.masked.fasta';genome.parent.mkdir(parents=True,exist_ok=True);genome.write_text('>chr1\n'+dna+'\n')
        reference(genome,aid,0,accepted/f'filtered_assemblies/{aid}')
        identity=json.loads((accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json').read_text())
        folder=accepted/f'structural_annotation/final/{aid}';folder.mkdir(parents=True);(folder/f'{aid}.gff3').write_text(model)
        readsets=[dict(readset_id=rid,library_id=rid,layout='single',strandedness='reverse',
            bam=f'readsets/{aid}/{rid}/aligned.bam',bam_index=f'readsets/{aid}/{rid}/aligned.bam.csi') for rid in ('R1','R2')] if aid=='A' else []
        row=dict(assembly_id=aid,mode=mode,status='ready',readsets=readsets)
        dump(folder/'annotation_status.json',dict(row,stage='final_annotation',reference_identity=identity['reference_identity'],gff3_sha256=checksum(folder/f'{aid}.gff3')))
        rows.append(row)
    dump(accepted/'structural_annotation/prediction_status.json',dict(status='ready',assemblies=rows))
    sam=root/'sams/A.sam';sam.parent.mkdir()
    with sam.open('w') as h:
        h.write('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:2000\n')
        for j,(start,end,strand) in enumerate([(200,300,'+'),(400,500,'+'),(1100,1200,'-')]):
            for offset in range(0,61,3):
                n=80-offset;pos=start-n;seq=dna[pos:start]+dna[end:end+100-n]
                for repeat in range(4):
                    h.write(f'j{j}_{offset}_{repeat}\t{16 if strand=="+" else 0}\tchr1\t{pos+1}\t60\t{n}M100N{100-n}M\t*\t0\t0\t{seq}\t'+ 'I'*100+f'\tNH:i:1\tNM:i:0\tXS:A:{strand}\n')


def run_case(root,name,resume=False,enabled=True):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'junction_qc.nf'),'-profile','slurm','-w',str(case/'work'),
        '--junction_inputs_root',str(root/'accepted'),'--outdir',str(case/'results'),
        '--portcullis_cpus','2','--portcullis_memory','8 GB','--run_junction_qc',str(enabled).lower()]
    if resume:command.append('-resume')
    print(f'Running actual Portcullis ({name})...',flush=True);execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    native=[t for t in tasks if 'RUN_PORTCULLIS (' in t['name']]
    if not enabled:
        assert not native;print('PASS: disabled',flush=True);return dict(case=name,passed=True)
    assert len(native)==2
    if name in ('resume','report_change','model_change'):assert all(t['status']=='CACHED' for t in native),native
    for aid in ('A','B'):
        out=case/f'results/qc/junctions/{aid}';state=json.loads((out/'junction_status.json').read_text())
        assert state['status']=='ready' and state['assessed']==(aid=='A')
        assert state['model_sha256']==checksum(root/f'accepted/structural_annotation/final/{aid}/{aid}.gff3')
        with (out/'transcript_junction_support.tsv').open() as h:rows={r['transcript_id']:r for r in csv.DictReader(h,delimiter='\t')}
        assert rows['t3']['all_pass_fraction']=='.'
        if aid=='A':
            with (out/'junction_support.tsv').open() as h:junctions=list(csv.DictReader(h,delimiter='\t'))
            assert {(r['intron_start'],r['intron_end'],r['strand']) for r in junctions}=={('201','300','+'),('401','500','+'),('1101','1200','-')}
            assert all(r['status']=='pass' for r in junctions),junctions
            assert all(r['passing_readsets']=='R1,R2' for r in junctions)
        else:assert rows['t1']['status']=='no_rna' and rows['t1']['all_pass_fraction']=='.'
    print(f'PASS: {name}',flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual tiny Portcullis including native low-junction fallback, both strands, per-readset support, no-RNA/single-exon, cache and model-only changes; no mapping, BRAKER or repeats; full random-forest pilot pending')
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
        for name in ('main.nf','junction_qc.nf','nextflow.config'):shutil.copy2(ROOT/name,pipeline/name)
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        (root/'main_syntax.console.log').write_text(probe.stdout+probe.stderr)
        if probe.returncode==0 or 'Provide --sample_sheet' not in probe.stdout+probe.stderr:
            raise RuntimeError('main.nf syntax/input guard failed; inspect '+str(root/'main_syntax.console.log')+
                               '\n'+(probe.stdout+probe.stderr)[-6000:])
        summary['checks'].append(dict(case='main_syntax',passed=True))
        execute(['nextflow','-C',str(ROOT/'nextflow.config'),'run',str(ROOT/'tests/utr_fixture_bams.nf'),'-profile','slurm',
            '--fixture_sams',str(root/'sams'),'--outdir',str(root/'accepted'),'-w',str(root/'bam_work'),'-resume'],root,root/'bam_fixture.console.log')
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        with (pipeline/'py_scripts/junction_support.py').open('a') as h:h.write('\n# Report-only cache regression.\n')
        summary['checks'].append(run_case(root,'report_change',True))
        for aid in ('A','B'):
            folder=root/f'accepted/structural_annotation/final/{aid}';gff=folder/f'{aid}.gff3'
            gff.write_text(gff.read_text()+'# Model export change, no BAM change.\n')
            status=json.loads((folder/'annotation_status.json').read_text());status['gff3_sha256']=checksum(gff);dump(folder/'annotation_status.json',status)
        summary['checks'].append(run_case(root,'model_change',True))
        summary['checks'].append(run_case(root,'disabled',True,False));summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root / "acceptance.json"}',flush=True)
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
