#!/usr/bin/env python3
"""Native InGenAnnot on tiny prealigned evidence, without prediction or remapping."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from validate_evidence_aed_cluster import fixtures as aed_fixtures
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from evidence_aed import normalize
from evidence_aed_selection import score
from prepare_protein_evidence import checksum


def fixtures(root):
    aed_fixtures(root)
    accepted=root/'accepted'
    predictions=json.loads((accepted/'structural_annotation/prediction_status.json').read_text())
    for row in predictions['assemblies']:
        aid=row['assembly_id']; mode=row['mode']
        gff=accepted/f'structural_annotation/braker/{aid}/{aid}.gff3'
        text=gff.read_text().replace('exon\t501\t1403\t.\t+\t.\tParent=t1','exon\t501\t1403\t.\t+\t.\tID=e1;Parent=t1')
        text=text.replace('CDS\t501\t1403\t.\t+\t0\tParent=t1','CDS\t501\t1403\t.\t+\t0\tID=c1;Parent=t1')
        gff.write_text(text)
        row['gff3_sha256']=checksum(gff)
        validation=accepted/f'structural_annotation/braker_status/{aid}.prediction.json'
        validation.write_text(json.dumps(row))
        rna=gff.parent/'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff' if mode in ('rna','rna_protein') else None
        protein=None
        output=accepted/f'structural_annotation/evidence_aed/{aid}'
        output.mkdir(parents=True)
        if mode in ('protein','rna_protein'):
            protein=output/f'protein_alignment/{aid}.miniprot.gff'
            protein.parent.mkdir()
            protein.write_text('##gff-version 3\nchr1\tminiprot\tmRNA\t501\t1400\t100\t+\t.\tID=MP1;Target=reference_protein 1 300;Identity=1.0\n'
                               'chr1\tminiprot\tCDS\t501\t1400\t100\t+\t0\tParent=MP1;Target=reference_protein 1 300\n')
        identity=json.loads((accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json').read_text())
        normalize(identity,rna,protein,root/f'normalized_{aid}')
        score(gff,identity,root/f'normalized_{aid}',output,row)
    (accepted/'structural_annotation/prediction_status.json').write_text(json.dumps(predictions))


def run_case(root,suffix,resume=False,penalty=.25,blocked=False):
    case=root/'pilot'; case.mkdir(exist_ok=True)
    pipeline=root/'pipeline'
    command=['nextflow','-log',str(case/f'{suffix}.nextflow.log'),'-C',str(pipeline/'nextflow.config'),
             'run',str(pipeline/'ingenannot_benchmark.nf'),'-profile','slurm','-w',str(case/'work'),
             '--aed_inputs_root',str(root/'accepted'),'--outdir',str(case/'results'),
             '--ingenannot_rna_penalty',str(penalty),'--ingenannot_aed_cpus','2',
             '--ingenannot_aed_memory','4 GB']
    if resume: command.append('-resume')
    (case/f'{suffix}.command.json').write_text(json.dumps(command,indent=2))
    print(f'Running native InGenAnnot ({suffix})...',flush=True)
    with (case/f'{suffix}.console.log').open('w') as log:
        result=subprocess.run(command,cwd=case,stdout=log,stderr=subprocess.STDOUT)
    console=(case/f'{suffix}.console.log').read_text(errors='replace')
    if blocked:
        if not result.returncode or 'evidence changed since current AED' not in console:
            raise RuntimeError('Changed evidence was not blocked as expected\n'+console[-6000:])
    else:
        if result.returncode: raise RuntimeError(f'{suffix}: exit {result.returncode}\n'+console[-6000:])
        for aid in 'ABCD':
            folder=case/f'results/structural_annotation/ingenannot_benchmark/{aid}'
            status=json.loads((folder/'benchmark_status.json').read_text())
            assert status['status']=='ready' and status['models_unchanged'] and status['transcripts']==1
            assert status['profile']['rna_overflow_penalty']==penalty
            with (folder/'benchmark_scores.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
            assert row['transcript_id']=='t1' and not any('combined' in k for k in row)
            assert (row['ingenannot_rna_aed']!='.')==(aid in 'AB')
            assert (row['ingenannot_protein_aed']!='.')==(aid in 'AC')
            if aid in 'AB': assert float(row['ingenannot_rna_aed'])==0
            if aid in 'AC': assert 0<=float(row['ingenannot_protein_aed'])<.05
        with (case/'results/pipeline/pipeline_trace.txt').open() as handle: tasks=list(csv.DictReader(handle,delimiter='\t'))
        assert not any(re.search(r'(ALIGN_AED_PROTEINS|BRAKER|REPEAT_MODEL|ALIGN_RNA_READSET)',t['name']) for t in tasks)
        for name in ['PREPARE_INGENANNOT_EVIDENCE','RUN_INGENANNOT_AED','REPORT_INGENANNOT_BENCHMARK']:
            found=[t for t in tasks if t['name'].startswith(name+' (')]
            assert len(found)==4
            if suffix=='resume' or (suffix=='profile_change' and name=='PREPARE_INGENANNOT_EVIDENCE') or (suffix=='report_change' and name!='REPORT_INGENANNOT_BENCHMARK'):
                assert all(t['status']=='CACHED' for t in found),f'{name} not cached: {found}'
            if suffix=='report_change' and name=='REPORT_INGENANNOT_BENCHMARK':
                assert all(t['status']=='COMPLETED' for t in found),f'Report edit did not rerun reporting: {found}'
    print(f'PASS: native InGenAnnot {suffix}',flush=True)
    return dict(case=suffix,passed=True,exit_code=result.returncode)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--outdir',required=True); p.add_argument('--resume',action='store_true'); a=p.parse_args()
    if sys.flags.optimize: p.error('Run without -O/PYTHONOPTIMIZE')
    root=Path(a.outdir).resolve(); root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()) and not a.resume: p.error('--outdir must be empty')
    if a.resume and not (root/'accepted/preflight/resolved_inputs.json').is_file(): p.error('--resume requires original fixtures')
    summary=dict(status='failed',checks=[],scope='Native InGenAnnot AED on tiny prealigned evidence, four modes, matched baseline, resume/profile-change cache and changed-source blocking; UTR refinement not run')
    try:
        version=subprocess.run(['nextflow','-version'],capture_output=True,text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b',version.stdout+version.stderr): raise RuntimeError('Load Nextflow 23.10.1')
        with (root/'unit_tests.log').open('w') as log:
            tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-p','test_*.py'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if tests.returncode: raise RuntimeError('Python tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests',passed=True))
        if not a.resume: fixtures(root)
        # A private source copy lets this test edit the reporting helper without
        # touching the user's checkout or a concurrently running workflow.
        pipeline=root/'pipeline'
        for relative in ['nextflow.config','ingenannot_benchmark.nf','modules/ingenannotBenchmark.nf',
                         'environments/ingenannot.yaml','py_scripts/ingenannot_inputs.py',
                         'py_scripts/ingenannot_benchmark.py','py_scripts/evidence_aed.py',
                         'py_scripts/prepare_protein_evidence.py']:
            destination=pipeline/relative
            destination.parent.mkdir(parents=True,exist_ok=True)
            if not destination.exists(): shutil.copy2(ROOT/relative,destination)
        summary['checks'].append(run_case(root,'recover' if a.resume else 'fresh',a.resume))
        summary['checks'].append(run_case(root,'resume',True))
        report_helper=pipeline/'py_scripts/ingenannot_benchmark.py'
        with report_helper.open('a') as handle: handle.write('\n# Report-only cache regression fixture.\n')
        summary['checks'].append(run_case(root,'report_change',True))
        summary['checks'].append(run_case(root,'profile_change',True,penalty=0))
        source=root/'accepted/structural_annotation/evidence_aed/A/protein_alignment/A.miniprot.gff'
        original=source.read_bytes()
        try:
            source.write_bytes(original+b'# changed\n')
            summary['checks'].append(run_case(root,'changed_evidence',True,blocked=True))
        finally: source.write_bytes(original)
        summary['status']='passed'
    except (AssertionError,ValueError,OSError,RuntimeError) as exc:
        summary['error']=str(exc); print(f'FAIL: {exc}',flush=True)
    (root/'acceptance.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'Acceptance summary: {root/"acceptance.json"}')
    return 0 if summary['status']=='passed' else 1


if __name__=='__main__': raise SystemExit(main())
