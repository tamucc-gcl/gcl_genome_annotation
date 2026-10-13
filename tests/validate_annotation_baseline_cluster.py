#!/usr/bin/env python3
"""Actual AGAT baseline contracts on tiny predictions; no BRAKER or database downloads.

BUSCO and OMArk correctness/runtime acceptance is a separate real-data pilot.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

from validate_preflight_cluster import fixtures as preflight_fixtures, table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'py_scripts'))
from prepare_protein_evidence import checksum
from reference_repeats import fasta


def fixtures(root):
    folder = root / 'fixtures'
    preflight_fixtures(folder)
    (folder / 'genome.fa').write_text('>s\n' + 'ATG' * 100 + '\n')
    header = ['assembly_id', 'assembly_fasta', 'use_protein_evidence', 'run_busco', 'run_omark']
    rows = [[aid, 'genome.fa', 'false', 'false', 'false'] for aid in ('A', 'B', 'C')]
    table(folder / 'baseline.csv', header, rows[:2])
    table(folder / 'added.csv', header, rows)
    accepted = folder / 'accepted'
    (accepted / 'repeats').mkdir(parents=True)
    validation_dir = accepted / 'structural_annotation/braker_status'
    validation_dir.mkdir(parents=True)
    status = []
    for aid in ('A', 'B', 'C'):
        (accepted / 'repeats' / f'{aid}.masked.fasta').write_text((folder / 'genome.fa').read_text())
        prediction = accepted / 'structural_annotation/braker' / aid
        prediction.mkdir(parents=True)
        gff, proteins = prediction / f'{aid}.gff3', prediction / f'{aid}_proteins.fa'
        lines = ['##gff-version 3']
        for gid, start, end, tids in [('g1', 1, 90, ['t1', 't2']), ('g2', 121, 210, ['t3'])]:
            lines.append(f's\tBRAKER\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}')
            for tid in tids:
                lines += [f's\tBRAKER\tmRNA\t{start}\t{end}\t.\t+\t.\tID={tid};Parent={gid}',
                          f's\tBRAKER\texon\t{start}\t{end}\t.\t+\t.\tID={tid}.exon;Parent={tid}',
                          f's\tBRAKER\tCDS\t{start}\t{end}\t.\t+\t0\tID={tid}.cds;Parent={tid}']
        gff.write_text('\n'.join(lines) + '\n')
        proteins.write_text(''.join(f'>{tid}\n' + 'M' * 30 + '\n' for tid in ('t1', 't2', 't3')))
        row = dict(assembly_id=aid, status='ready', mode='ab_initio', genes=2, coding_transcripts=3,
                   gff3_sha256=checksum(gff), proteins_sha256=checksum(proteins))
        contigs = [{'id': header.split()[0], 'length': len(sequence),
                    'uppercase_sha256': hashlib.sha256(sequence.upper().encode('ascii')).hexdigest()}
                   for header, sequence in fasta(folder / 'genome.fa')]
        row['reference_identity'] = hashlib.sha256(json.dumps(contigs, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        (validation_dir / f'{aid}.prediction.json').write_text(json.dumps(row))
        status.append(row)
    (accepted / 'structural_annotation/prediction_status.json').write_text(json.dumps(dict(status='ready', assemblies=status)))


def run_case(root, suffix, sheet='baseline.csv', resume=False, blocked=False):
    case = root / 'baseline'
    case.mkdir(exist_ok=True)
    folder = root / 'fixtures'
    output = case / 'results'
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / 'annotation_qc.nf'), '-profile', 'slurm', '-w', str(case / 'work'),
               '--sample_sheet', str(folder / sheet), '--taxid', '999', '--taxdump_dir', str(folder / 'taxdump'),
               '--baseline_inputs_root', str(folder / 'accepted'), '--outdir', str(output)]
    if resume:
        command.append('-resume')
    (case / f'{suffix}.command.json').write_text(json.dumps(command, indent=2))
    print(f'Running annotation baseline ({suffix}; actual AGAT)...', flush=True)
    with (case / f'{suffix}.console.log').open('w') as log:
        result = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    if blocked:
        if result.returncode == 0:
            raise RuntimeError('Changed prediction checksum did not block QC')
        if 'Prediction files changed since validation' not in (case / f'{suffix}.console.log').read_text(errors='replace'):
            raise RuntimeError('Blocked case failed for an unexpected reason; inspect its console log')
        print('PASS: changed prediction blocked', flush=True)
        return dict(case=suffix, passed=True, exit_code=result.returncode)
    status = output / 'qc/baseline/baseline_status.json'
    if result.returncode or not status.is_file():
        tail = '\n'.join((case / f'{suffix}.console.log').read_text(errors='replace').splitlines()[-50:])
        raise RuntimeError(f'{suffix}: exit {result.returncode}; missing/failed baseline.\n{tail}')
    rows = json.loads(status.read_text())['assemblies']
    expected = {'A', 'B', 'C'} if sheet == 'added.csv' else {'A', 'B'}
    assert {r['assembly_id'] for r in rows} == expected
    for row in rows:
        assert row['genes'] == 2 and row['coding_transcripts'] == 3 and row['representative_proteins'] == 2
        assert row['analyses'] == dict(agat='ready', busco='not_requested', omark='not_requested')
        aid = row['assembly_id']
        assert (output / f'qc/baseline/{aid}/selection/representatives.tsv').is_file()
        assert (output / f'qc/baseline/{aid}/agat/{aid}_annotation_stats.yaml').stat().st_size > 0
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    assert not any(re.search(r':(BUSCO|OMARK|DOWNLOAD_BUSCO_DB|DOWNLOAD_OMARK_DB)( |$)', r['name']) for r in trace)
    assert not any('BRAKER' in r['name'] for r in trace)
    if suffix in ('resume', 'added'):
        tasks = [r for r in trace if re.search(r':(PREPARE_ANNOTATION_BASELINE|AGAT) \((A|B)\)', r['name'])]
        assert len(tasks) == 4 and all(r['status'] == 'CACHED' for r in tasks), f'Unexpected rerun: {tasks}'
    print(f'PASS: baseline {suffix}', flush=True)
    return dict(case=suffix, passed=True, exit_code=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O/PYTHONOPTIMIZE')
    root = Path(args.outdir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()) and not args.resume:
        parser.error('--outdir must be empty')
    if args.resume and not (root / 'fixtures/baseline.csv').is_file():
        parser.error('--resume requires the original fixtures')
    summary = dict(status='failed', checks=[], scope='Actual AGAT, multi-assembly baseline, representative selection, disabled metrics, resume, additions and changed-file blocking; BUSCO/OMArk real pilot pending')
    try:
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version.stdout + version.stderr):
            raise RuntimeError('Load Nextflow 23.10.1')
        with (root / 'unit_tests.log').open('w') as log:
            tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'), '-p', 'test_*.py'],
                                   cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if tests.returncode:
            raise RuntimeError('Python tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests', passed=True))
        if not args.resume:
            fixtures(root)
        summary['checks'].append(run_case(root, 'recover' if args.resume else 'fresh', resume=args.resume))
        summary['checks'].append(run_case(root, 'resume', resume=True))
        summary['checks'].append(run_case(root, 'added', 'added.csv', True))
        path = root / 'fixtures/accepted/structural_annotation/braker_status/A.prediction.json'
        original = path.read_bytes()
        try:
            data = json.loads(original); data['proteins_sha256'] = 'wrong'
            path.write_text(json.dumps(data))
            summary['checks'].append(run_case(root, 'changed_prediction', resume=True, blocked=True))
        finally:
            path.write_bytes(original)
        summary['status'] = 'passed'
    except (AssertionError, ValueError, OSError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', flush=True)
    (root / 'acceptance.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {root / "acceptance.json"}')
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
