#!/usr/bin/env python3
"""Actual before/after AGAT and final-ID structural report; no prediction/searches."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from validate_selection_cluster import fixtures as selection_fixtures, MODES
from validate_utr_cluster import execute,dump,ROOT
sys.path.insert(0,str(ROOT/'py_scripts'))
from model_selection import select
from prepare_protein_evidence import checksum


def fixtures(root):
    selection_fixtures(root);accepted=root/'accepted';rows=[]
    for aid,mode in MODES.items():
        genome=accepted/f'repeats/{aid}.masked.fasta';genome.parent.mkdir(exist_ok=True)
        shutil.copy2(root/f'source_{aid}/genome.fa',genome)
        identity=accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json';ref=json.loads(identity.read_text())
        ref['reference_identity']=hashlib.sha256(json.dumps(ref['contigs'],sort_keys=True,separators=(',',':')).encode()).hexdigest();dump(identity,ref)
        folder=accepted/f'structural_annotation/final/{aid}';v=folder/'annotation_status.json';state=json.loads(v.read_text())
        state.update(reference_identity=ref['reference_identity'],genes=2,coding_transcripts=3);dump(v,state)
        assessment=accepted/f'structural_annotation/assessment/{aid}';ast=assessment/'assessment_status.json'
        state=json.loads(ast.read_text());state['reference_identity']=ref['reference_identity'];dump(ast,state)
        select(aid,folder/f'{aid}.gff3',folder/f'{aid}_proteins.fa',folder/f'{aid}_cds.fa',v,identity,
            assessment/'transcript_assessment.tsv',ast,accepted/f'structural_annotation/selection/{aid}')
        rows.append(dict(assembly_id=aid,settings=dict(run_busco=False,run_omark=False,busco_lineage=None,taxid=141886)))
    dump(accepted/'preflight/resolved_inputs.json',dict(status='ready',assemblies=rows))


def run_case(root,name,resume=False):
    case=root/'connected';case.mkdir(exist_ok=True);pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{name}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
        'run',str(pipeline/'structural_acceptance.nf'),'-profile','slurm','-w',str(case/'work'),
        '--structural_inputs_root',str(root/'accepted'),'--outdir',str(case/'results')]
    if resume:command.append('-resume')
    if name=='prefix_change':command+=['--final_id_prefix','gcl']
    if name=='disabled':command+=['--run_finalization','false']
    print(f'Running final structural QC ({name}; actual before/after AGAT)...',flush=True)
    execute(command,case,case/f'{name}.console.log')
    with (case/'results/pipeline/pipeline_trace.txt').open() as handle:tasks=list(csv.DictReader(handle,delimiter='\t'))
    if name=='disabled':
        assert not tasks,tasks
        print('PASS: disabled',flush=True);return dict(case=name,passed=True)
    assert len([t for t in tasks if 'RENAME_GENE_MODELS (' in t['name']])==4
    assert len([t for t in tasks if re.search(r':AGAT \(',t['name'])])==8
    for t in tasks:
        if name=='resume':assert t['status']=='CACHED',t
        if name=='report_change' and ('STRUCTURAL_QC (' not in t['name'] and 'STRUCTURAL_QC_INDEX' not in t['name']):
            assert t['status']=='CACHED',t
        if name=='prefix_change' and 'FINAL_ANNOTATION_QC:' in t['name']:assert t['status']=='CACHED',t
        assert not re.search(r'(REPEAT_MODEL|ALIGN_RNA_READSET|RUN_INGENANNOT_AED|:BRAKER \(|:BUSCO \(|:OMARK \()',t['name']),t
    for aid in MODES:
        final=case/f'results/structural_annotation/identified/{aid}'
        state=json.loads((final/'annotation_status.json').read_text());assert state['genes']==1 and state['coding_transcripts']==2
        prefix=f'gcl_{aid}' if name=='prefix_change' else aid
        with (final/'transcript_assessment.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        assert {r['transcript_id'] for r in rows}=={f'{prefix}_g1.t1',f'{prefix}_g1.t2'}
        with (final/'transcript_id_map.tsv').open() as handle:mapped=list(csv.DictReader(handle,delimiter='\t'))
        assert len(mapped)==3 and next(r for r in mapped if r['transcript_id']=='t3')['final_transcript_id']=='.'
        qc=case/f'results/qc/structural/{aid}';state=json.loads((qc/'structural_qc_status.json').read_text())
        assert [r['transcripts'] for r in state['model_sets']]==[3,2]
        assert state['functional_annotation']=='pending' and not state['confidence_list_created']
        assert len(list((qc/'plots').glob('*.svg')))==10 and len(list((qc/'plots').glob('*.png')))==10
        assert state['native_qc_sources']['before'] and state['native_qc_sources']['after']
        assert (qc/'README.md').is_file() and not list(qc.rglob('*.html'))
        for filename,sha in state['outputs'].items():assert checksum(qc/filename)==sha,filename
        selected=json.loads((case/f'results/qc/selected/{aid}/selection/selection.json').read_text())
        assert selected['model_set']=='selected' and selected['representative_proteins']==1
    assert (case/'results/qc/structural/README.md').is_file()
    print(f'PASS: {name}',flush=True);return dict(case=name,passed=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if sys.flags.optimize:p.error('Run without -O')
    root=Path(a.outdir).resolve();root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume:p.error('--outdir must be empty')
    summary=dict(status='failed',checks=[],scope='Actual Nextflow final IDs, all four synthetic accepted evidence modes, before/after AGAT, complete maps, selected-QC representative provenance, ten PNG/SVG report figures, resume/report-only cache and prefix change; no prediction or searches; full BUSCO/OMArk pilot pending')
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
        for name in ('main.nf','final_structural_qc.nf','nextflow.config'):shutil.copy2(ROOT/name,pipeline/name)
        shutil.copy2(ROOT/'tests/final_structural_fixture.nf',pipeline/'structural_acceptance.nf')
        probe=subprocess.run(['nextflow','-C',str(pipeline/'nextflow.config'),'run',str(pipeline/'main.nf')],cwd=pipeline,capture_output=True,text=True)
        (root/'main_syntax.console.log').write_text(probe.stdout+probe.stderr)
        if probe.returncode==0 or 'Provide --sample_sheet' not in probe.stdout+probe.stderr:raise RuntimeError('main syntax/input guard failed\n'+(probe.stdout+probe.stderr)[-6000:])
        summary['checks'].append(dict(case='main_syntax',passed=True))
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        with (pipeline/'py_scripts/structural_report.py').open('a') as handle:handle.write('\n# Report cache regression.\n')
        for name in ('report_change','prefix_change','disabled'):summary['checks'].append(run_case(root,name,True))
        summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc);print(f'FAIL: {exc}',flush=True)
    dump(root/'acceptance.json',summary);print(f'Acceptance summary: {root / "acceptance.json"}',flush=True)
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__':raise SystemExit(main())
