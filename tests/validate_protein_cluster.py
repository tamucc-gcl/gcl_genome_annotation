#!/usr/bin/env python3
"""Cluster acceptance: local protein fixtures, disabled mode, failures and resume."""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from validate_preflight_cluster import fixtures, table

ROOT = Path(__file__).resolve().parents[1]


def run_case(folder, fixture_dir, name, sheet, expected, source_count, resume=False, added=0):
    case = folder / name
    case.mkdir(exist_ok=True)
    output = case / 'results'
    output.mkdir(exist_ok=True)
    suffix = 'added_assembly' if added else ('resume' if resume else 'fresh')
    command = ['nextflow', '-log', str(case / (suffix + '.nextflow.log')),
        '-C', str(ROOT / 'nextflow.config'), 'run', str(ROOT / 'protein_evidence.nf'),
        '-profile', 'slurm', '-w', str(case / 'work'), '--sample_sheet', str(fixture_dir / sheet),
        '--taxid', '999', '--taxdump_dir', str(fixture_dir / 'taxdump'), '--outdir', str(output)]
    if resume:
        command.append('-resume')
    print(f'Running {name} ({suffix}; local synthetic proteins only)...', flush=True)
    with (case / (suffix + '.console.log')).open('w') as log:
        process = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    status_path = output / 'protein_evidence/protein_status.json'
    if not status_path.is_file():
        tail = '\n'.join((case / (suffix + '.console.log')).read_text(errors='replace').splitlines()[-35:])
        raise RuntimeError(f'{name}: exit {process.returncode}, missing report.\n{tail}')
    status = json.loads(status_path.read_text())
    states = {row['assembly_id']: row['status'] for row in status['assemblies']}
    assert states == expected, f'{name}: unexpected states {states}'
    success = all(value in ('ready', 'not_requested') for value in expected.values())
    assert (process.returncode == 0) == success, f'{name}: exit {process.returncode}'
    source_manifests = list((output / 'protein_sources').glob('*/source.json'))
    assert len(source_manifests) == source_count, f'{name}: wrong unique source count'
    for row in status['assemblies']:
        proteins = output / 'protein_evidence/assemblies' / row['assembly_id'] / 'proteins.fa'
        if row['status'] == 'ready':
            assert proteins.is_file() and row['protein_count'] > 0
        else:
            assert not proteins.exists(), f'{name}: absent/blocked evidence produced a FASTA'
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    sources = [row for row in trace if row['name'].startswith('ACQUIRE_PROTEIN_SOURCE')]
    assert len(sources) == source_count, f'{name}: acquisition task count differs from unique sources'
    if resume:
        assert all(row['status'] == 'CACHED' for row in sources), 'Unchanged source acquisitions were not cached'
        combines = [row for row in trace if row['name'].startswith('COMBINE_PROTEIN_EVIDENCE')]
        assert len(combines) == len(expected), 'Missing per-assembly combination task'
        assert sum(row['status'] == 'CACHED' for row in combines) == len(expected) - added, (
            'Unchanged per-assembly combinations were not cached')
        assert sum(row['status'] == 'COMPLETED' for row in combines) == added, (
            'Expected new assembly combinations were not run')
    print(f'PASS: {name} ({suffix})', flush=True)
    return {'case': name, 'resume': resume, 'passed': True, 'exit_code': process.returncode,
            'assembly_states': states, 'source_count': source_count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, help='New empty shared cluster directory')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O/PYTHONOPTIMIZE')
    folder = Path(args.outdir).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if any(folder.iterdir()):
        parser.error('--outdir must be empty')
    summary = {'status': 'failed', 'checks': []}
    try:
        with (folder / 'unit_tests.log').open('w') as log:
            tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'),
                '-p', 'test_*.py', '-v'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if tests.returncode:
            raise RuntimeError(f'Protein unit tests failed; inspect {folder / "unit_tests.log"}')
        summary['checks'].append({'case': 'unit_tests', 'passed': True})
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version.stdout + version.stderr):
            raise RuntimeError('Load Nextflow 23.10.1')
        fixture_dir = folder / 'fixtures'
        fixtures(fixture_dir)
        (fixture_dir / 'one.fa').write_text('>shared from source one\nMPEPTIDE*\n')
        (fixture_dir / 'two.fa').write_text('>shared from source two\nMKL\n')
        (fixture_dir / 'bad.fa').write_text('>internal_stop\nMP*EP\n')
        header = ['assembly_id', 'assembly_fasta', 'use_protein_evidence', 'orthodb_clade', 'protein_fastas']
        table(fixture_dir / 'ready_proteins.csv', header, [
            ['A', 'genome.fa', 'true', 'none', 'one.fa;two.fa'],
            ['B', 'genome.fa', 'true', 'none', 'one.fa'],
            ['disabled', 'genome.fa', 'false', 'none', 'none']])
        table(fixture_dir / 'blocked_proteins.csv', header, [
            ['bad_a', 'genome.fa', 'true', 'none', 'bad.fa'],
            ['bad_b', 'genome.fa', 'true', 'none', 'bad.fa'],
            ['independent', 'genome.fa', 'true', 'none', 'one.fa']])
        table(fixture_dir / 'disabled_proteins.csv', header,
              [['disabled', 'genome.fa', 'false', 'none', 'none']])
        table(fixture_dir / 'added_proteins.csv', header, [
            ['A', 'genome.fa', 'true', 'none', 'one.fa;two.fa'],
            ['B', 'genome.fa', 'true', 'none', 'one.fa'],
            ['disabled', 'genome.fa', 'false', 'none', 'none'],
            ['C', 'genome.fa', 'true', 'none', 'two.fa']])
        ready = {'A': 'ready', 'B': 'ready', 'disabled': 'not_requested'}
        summary['checks'].append(run_case(folder, fixture_dir, 'ready', 'ready_proteins.csv', ready, 2))
        summary['checks'].append(run_case(folder, fixture_dir, 'ready', 'ready_proteins.csv', ready, 2, resume=True))
        summary['checks'].append(run_case(folder, fixture_dir, 'ready', 'added_proteins.csv',
            {**ready, 'C': 'ready'}, 2, resume=True, added=1))
        summary['checks'].append(run_case(folder, fixture_dir, 'blocked', 'blocked_proteins.csv',
            {'bad_a': 'blocked', 'bad_b': 'blocked', 'independent': 'ready'}, 2))
        summary['checks'].append(run_case(folder, fixture_dir, 'disabled', 'disabled_proteins.csv',
            {'disabled': 'not_requested'}, 0))
        summary['status'] = 'passed'
    except (OSError, ValueError, AssertionError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', file=sys.stderr)
    (folder / 'acceptance.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {folder / "acceptance.json"}', flush=True)
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
