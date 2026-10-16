#!/usr/bin/env python3
"""Tiny native UTR ranking/refinement, gffread and AGAT acceptance; no prediction."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from prepare_protein_evidence import checksum


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value,indent=2)+'\n')


def fixtures(root):
    accepted=root/'accepted'; predictions=[]
    dna='A'*2000
    contigs=[dict(id='chr1',length=len(dna),uppercase_sha256=hashlib.sha256(dna.encode()).hexdigest())]
    identity=dict(reference_identity=hashlib.sha256(json.dumps(contigs,sort_keys=True,separators=(',',':')).encode()).hexdigest(),contigs=contigs)
    for aid,strand,mode in [('A','+','rna_protein'),('B','+','protein'),('C','-','rna')]:
        gff=accepted/f'structural_annotation/braker/{aid}/{aid}.gff3'; gff.parent.mkdir(parents=True)
        gff.write_text('##gff-version 3\n'+''.join(f'chr1\tBRAKER\t{k}\t{s}\t{e}\t.\t{strand}\t{phase}\t{attrs}\n'
            for k,s,e,phase,attrs in [('gene',501,900,'.','ID=g1'),('mRNA',501,900,'.','ID=t1;Parent=g1'),
                ('exon',501,650,'.','ID=e1;Parent=t1'),('exon',751,900,'.','ID=e2;Parent=t1'),
                ('CDS',501,650,'0','ID=c1;Parent=t1'),('CDS',751,900,'0','ID=c2;Parent=t1')]))
        genome=accepted/f'repeats/{aid}.masked.fasta';genome.parent.mkdir(exist_ok=True);genome.write_text('>chr1\n'+dna+'\n')
        dump(accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json',identity)
        readsets=[]
        if mode!='protein':
            readsets=[dict(readset_id='R1',layout='paired',strandedness='reverse',bam=f'readsets/{aid}/R1/aligned.bam',bam_index=f'readsets/{aid}/R1/aligned.bam.csi')]
            rna=gff.parent/'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff';rna.parent.mkdir(parents=True)
            rna.write_text(''.join(f'chr1\tStringTie\t{k}\t{s}\t{e}\t.\t{strand}\t.\tgene_id "r1"; transcript_id "rna1";\n'
                for k,s,e in [('transcript',451,950),('exon',451,650),('exon',751,950)]))
            sam=root/f'sams/{aid}.sam';sam.parent.mkdir(exist_ok=True)
            def alignment(offset):
                pos=451+offset if offset<200 else 751+offset-200
                cigar=f'{200-offset}M100N{offset+100-200}M' if offset<200<offset+100 else '100M'
                return pos,cigar
            with sam.open('w') as h:
                h.write('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:2000\n')
                for i,offset in enumerate(range(0,201,5)):
                    for mate,delta,flag in [(1,100,83 if strand=='+' else 99),(2,0,163 if strand=='+' else 147)]:
                        pos,cigar=alignment(offset+delta);other,_=alignment(offset+(0 if mate==1 else 100))
                        h.write(f'p{i}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t=\t{other}\t0\t'+ 'A'*100+'\t'+'I'*100+'\n')
        row=dict(assembly_id=aid,status='ready',mode=mode,readsets=readsets,gff3_sha256=checksum(gff),reference_identity=identity['reference_identity'])
        predictions.append(row); dump(accepted/f'structural_annotation/braker_status/{aid}.prediction.json',row)
    dump(accepted/'structural_annotation/prediction_status.json',dict(status='ready',assemblies=predictions))


def execute(command,cwd,log):
    with log.open('w') as h: result=subprocess.run(command,cwd=cwd,stdout=h,stderr=subprocess.STDOUT)
    if result.returncode: raise RuntimeError(f'{log}: exit {result.returncode}\n'+log.read_text(errors='replace')[-6000:])


def run_case(root,name,resume=False,enabled=True):
    case=root/'pilot';case.mkdir(exist_ok=True); pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'utr_refinement.nf'),
        '-profile','slurm','-w',str(case/'work'),'--utr_inputs_root',str(root/'accepted'),'--outdir',str(case/'results'),
        '--utr_ranking_cpus','2','--utr_ranking_memory','4 GB','--utr_refine',str(enabled).lower()]
    if resume: command.append('-resume')
    print(f'Running native UTR pilot ({name})...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    for aid in ('A','B','C'):
        folder=case/f'results/structural_annotation/utr/{aid}'
        status=json.loads((folder/'utr_status.json').read_text())
        assert status['status']=='ready' and status['protein_sequences_unchanged'] and status['cds_coordinates_and_phases_unchanged']
        assert status['changed_transcripts']==(1 if enabled and aid!='B' else 0)
        assert status['skipped']==(not enabled or aid=='B')
        with (folder/'utr_changes.tsv').open() as h: row=next(csv.DictReader(h,delimiter='\t'))
        if enabled and aid!='B': assert row['five_prime_utr_bp']=='50' and row['three_prime_utr_bp']=='50'
        assert (case/f'results/qc/utr/{aid}/agat/{aid}_annotation_stats').is_file()
    with (case/'results/pipeline/pipeline_trace.txt').open() as h: tasks=list(csv.DictReader(h,delimiter='\t'))
    ranked=[t for t in tasks if t['name'].startswith('RUN_UTR_REFINEMENT (')]
    assert len(ranked)==3
    if name in ('resume','report_change'): assert all(t['status']=='CACHED' for t in ranked),f'Ranking not cached: {ranked}'
    assert not any(re.search(r'(REPEAT_MODEL|ALIGN_RNA_READSET|BRAKER\s*\()',t['name']) for t in tasks)
    print(f'PASS: native UTR {name}',flush=True)
    return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize: p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Tiny actual InGenAnnot cluster/rank/refine, both strands, no-RNA/disabled pass-through, gffread invariance, AGAT, ranking cache; full-genome pilot and final AED rescoring pending')
    try:
        version=subprocess.run(['nextflow','-version'],capture_output=True,text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b',version.stdout+version.stderr):raise RuntimeError('Load Nextflow 23.10.1')
        with (root/'unit_tests.log').open('w') as log:
            tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if tests.returncode: raise RuntimeError('Python tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume: fixtures(root)
        files=['nextflow.config','utr_refinement.nf','modules/utrRefinement.nf','modules/evidenceAed.nf','modules/agat.nf',
            'environments/ingenannot.yaml','environments/agat.yaml','py_scripts/utr_refinement.py','py_scripts/evidence_aed.py',
            'py_scripts/prepare_protein_evidence.py','py_scripts/reference_repeats.py']
        for relative in files:
            dest=root/'pipeline'/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            if not dest.exists() or dest.read_bytes()!=(ROOT/relative).read_bytes():shutil.copy2(ROOT/relative,dest)
        execute(['nextflow','-C',str(ROOT/'nextflow.config'),'run',str(ROOT/'tests/utr_fixture_bams.nf'),'-profile','slurm',
            '--fixture_sams',str(root/'sams'),'--outdir',str(root/'accepted'),'-w',str(root/'bam_work'),
            '--rna_alignment_cpus','2','--rna_alignment_memory','4 GB','-resume'],root,root/'bam_fixture.console.log')
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        with (root/'pipeline/py_scripts/utr_refinement.py').open('a') as h:h.write('\n# Reporting helper cache regression.\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'disabled',True,False))
        summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root/"acceptance.json"}')
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
