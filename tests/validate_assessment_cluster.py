#!/usr/bin/env python3
"""Actual SeqKit/DIAMOND/BEDTools model assessment; synthetic accepted evidence handoffs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys

from test_model_assessment import fixtures as model_fixture
from validate_utr_cluster import execute, dump, ROOT
sys.path.insert(0,str(ROOT/'py_scripts'))
from prepare_protein_evidence import checksum


def fixtures(root):
    source=root/'source';source.mkdir();g,p,c,genome,identity,validation,metrics,native,coverage=model_fixture(source)
    rng=random.Random(19)
    codons=dict(A='GCT',C='TGT',D='GAT',E='GAA',F='TTT',G='GGT',H='CAT',I='ATT',K='AAA',L='CTG',M='ATG',N='AAT',P='CCT',Q='CAA',R='CGT',S='TCT',T='ACT',V='GTT',W='TGG',Y='TAT')
    amino='M'+''.join(rng.choice(list(codons)) for _ in range(32));dna=''.join(codons[a] for a in amino)+'A'
    other='M'+''.join(rng.choice(list(codons)) for _ in range(32))
    p.write_text(f'>t1\n{amino}\n>t2\n{amino}\n>t3\n{other}\n')
    c.write_text(f'>t1\n{dna}\n>t2\n{dna}\n>t3\n'+ 'A'*100+'\n')
    accepted=root/'accepted';rows=[]
    for aid,mode in [('A','rna_protein'),('B','protein'),('C','rna'),('D','ab_initio')]:
        folder=accepted/f'structural_annotation/final/{aid}';folder.mkdir(parents=True)
        for src,name in [(g,f'{aid}.gff3'),(p,f'{aid}_proteins.fa'),(c,f'{aid}_cds.fa')]:shutil.copy2(src,folder/name)
        ref=accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json';ref.parent.mkdir(parents=True);shutil.copy2(identity,ref)
        target=accepted/f'repeats/{aid}.masked.fasta';target.parent.mkdir(exist_ok=True);shutil.copy2(genome,target)
        (target.parent/f'{aid}.repeat.gff3').write_text('chr1\tRepeatMasker\trepeat_region\t101\t175\t.\t.\t.\tID=r1\nchr1\tRepeatMasker\trepeat_region\t326\t400\t.\t.\t.\tID=r2\nchr1\tRepeatMasker\trepeat_region\t101\t175\t.\t.\t.\tID=r3\n')
        library=accepted/f'repeats/denovo/{aid}_curatedRepeats.fasta';library.parent.mkdir(exist_ok=True);library.write_text('>rpt1#LINE/L1\n'+dna+'\n')
        state=json.loads(validation.read_text());state.update(assembly_id=aid,mode=mode,proteins_sha256=checksum(p),cds_sha256=checksum(c))
        dump(folder/'annotation_status.json',state)
        rows.append(dict(assembly_id=aid,mode=mode,status='ready',readsets=[]))
        has_rna=mode in ('rna_protein','rna');has_protein=mode in ('rna_protein','protein')
        aed=accepted/f'structural_annotation/aed/{aid}';aed.mkdir(parents=True)
        (aed/'aed_scores.tsv').write_text('gene_id\ttranscript_id\trna_aed\tprotein_aed\trna_evidence_id\tprotein_evidence_id\trna_penalty\n'+''.join(
            f'{gid}\t{tid}\t{0.25 if has_rna else "."}\t{0.1 if has_protein else "."}\trna1\tprotein1\tno\n' for gid,tid in [('g1','t1'),('g1','t2'),('g2','t3')]))
        dump(aed/'aed_status.json',dict(status='ready',stage='native_aed',model_sha256=checksum(g),reference_identity='ref',
            evidence=dict(rna=dict(available=has_rna),protein=dict(available=has_protein))))
        junction=accepted/f'qc/junctions/{aid}';junction.mkdir(parents=True)
        table=[dict(gene_id=gid,transcript_id=tid,status='assessed' if has_rna else 'no_rna',all_junctions=1,
            all_pass=1 if has_rna else '.',all_fail=0 if has_rna else '.',all_not_observed=0 if has_rna else '.',all_unresolved_strand=0 if has_rna else '.',all_pass_fraction=1 if has_rna else '.',
            coding_junctions=1,coding_pass=1 if has_rna else '.',coding_fail=0 if has_rna else '.',coding_not_observed=0 if has_rna else '.',coding_unresolved_strand=0 if has_rna else '.',coding_pass_fraction=1 if has_rna else '.',
            utr_junctions=0,utr_pass=0 if has_rna else '.',utr_fail=0 if has_rna else '.',utr_not_observed=0 if has_rna else '.',utr_unresolved_strand=0 if has_rna else '.',utr_pass_fraction='.') for gid,tid in [('g1','t1'),('g1','t2'),('g2','t3')]]
        with (junction/'transcript_junction_support.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(table[0]),delimiter='\t');w.writeheader();w.writerows(table)
        dump(junction/'junction_status.json',dict(status='ready',stage='junction_qc',model_sha256=checksum(g),reference_identity='ref',assessed=has_rna))
    dump(accepted/'structural_annotation/prediction_status.json',dict(status='ready',assemblies=rows))


def run_case(root,name,resume=False,te=True,overlap=.5,evalue=1e-10,evidence=True):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'model_assessment.nf'),'-profile','slurm','-w',str(case/'work'),
        '--assessment_inputs_root',str(root/'accepted'),'--outdir',str(case/'results'),
        '--run_te_flagging',str(te).lower(),'--te_flag_overlap_frac',str(overlap),'--te_flag_evalue',str(evalue),
        '--run_aed',str(evidence).lower(),'--run_junction_qc',str(evidence).lower()]
    if resume:command.append('-resume')
    print(f'Running model assessment ({name}; actual SeqKit/DIAMOND/BEDTools)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    for aid,mode in [('A','rna_protein'),('B','protein'),('C','rna'),('D','ab_initio')]:
        out=case/f'results/structural_annotation/assessment/{aid}'
        state=json.loads((out/'assessment_status.json').read_text())
        assert state['status']=='ready' and state['models_removed']==0 and state['transcripts']==3 and state['genes']==2
        assert state['model_sha256']==checksum(root/f'accepted/structural_annotation/final/{aid}/{aid}.gff3')
        with (out/'transcript_assessment.tsv').open() as h:rows={r['transcript_id']:r for r in csv.DictReader(h,delimiter='\t')}
        assert (rows['t1']['rna_aed']!='.')==(evidence and mode in ('rna_protein','rna'))
        assert (rows['t1']['protein_aed']!='.')==(evidence and mode in ('rna_protein','protein'))
        if not evidence:assert rows['t1']['junction_status']=='disabled'
        assert rows['t3']['donor_acceptor_summary']=='GTAG:1'
        if te:
            native=json.loads((out/'te_search/search/search.json').read_text())
            assert native['hits']==(0 if name=='no_hits' else 2),native
            assert state['te_flagged_genes']==(0 if name=='no_hits' else 1),state
            assert abs(float(rows['t1']['te_gene_repeat_fraction'])-.5)<1e-12
            with (out/'te/te_gene_report.tsv').open() as h:genes={r['gene_id']:r for r in csv.DictReader(h,delimiter='\t')}
            assert genes['g1']['cds_bp']=='100' # Shared CDS is not counted twice across isoforms.
        else:assert rows['t1']['te_gene_flag']=='.' and not state['te_enabled']
    with (case/'results/pipeline/pipeline_trace.txt').open() as h:tasks=list(csv.DictReader(h,delimiter='\t'))
    native=[t for t in tasks if re.search(r':(BUILD_TE_REPEAT_DB|SEARCH_TE_PROTEINS) \(',t['name'])]
    metric_tasks=[t for t in tasks if 'PREPARE_MODEL_METRICS (' in t['name']]
    assert len(metric_tasks)==4
    if name in ('resume','report_change','overlap_change','no_hits','disabled','evidence_disabled'):
        assert all(t['status']=='CACHED' for t in metric_tasks),metric_tasks
    if not te:assert not native
    else:
        assert len(native)==8,native
        if name in ('resume','report_change','overlap_change','evidence_disabled'):assert all(t['status']=='CACHED' for t in native),native
        if name=='no_hits':assert all(t['status']=='CACHED' for t in native if 'BUILD_TE_REPEAT_DB (' in t['name'])
    assert not any(re.search(r'(REPEAT_MODEL|ALIGN_RNA_READSET|RUN_INGENANNOT_AED|:BRAKER \()',t['name']) for t in tasks)
    print(f'PASS: {name}',flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual SeqKit/DIAMOND/BEDTools; both strands, gene-union overlap/shared isoforms, four synthetic accepted evidence modes, no-hit/disabled, cache and report/overlap-only reruns; no prediction, repeats or evidence recalculation')
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
        for name in ('main.nf','model_assessment.nf','nextflow.config'):shutil.copy2(ROOT/name,pipeline/name)
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        (root/'main_syntax.console.log').write_text(probe.stdout+probe.stderr)
        if probe.returncode==0 or 'Provide --sample_sheet' not in probe.stdout+probe.stderr:raise RuntimeError('main syntax/input guard failed\n'+(probe.stdout+probe.stderr)[-6000:])
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        with (pipeline/'py_scripts/model_assessment.py').open('a') as h:h.write('\n# Report cache regression.\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'overlap_change',True,overlap=.75))
        summary['checks'].append(run_case(root,'no_hits',True,overlap=.75,evalue=1e-100))
        summary['checks'].append(run_case(root,'disabled',True,te=False))
        summary['checks'].append(run_case(root,'evidence_disabled',True,evidence=False));summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root / "acceptance.json"}',flush=True)
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
