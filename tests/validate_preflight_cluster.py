#!/usr/bin/env python3
"""Run input-contract tests and tiny Nextflow preflights on the cluster."""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def table(path, header, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def fixtures(folder):
    folder.mkdir()
    (folder / 'genome.fa').write_text('>chr1\nACGTACGT\n')
    (folder / 'r1.fq').write_text('@read/1\nACGT\n+\nIIII\n')
    (folder / 'r2.fq').write_text('@read/2\nTGCA\n+\nIIII\n')
    tax = folder / 'taxdump'
    tax.mkdir()
    nodes = [('1', '1'), ('2759', '1'), ('7742', '2759'), ('7898', '7742'), ('999', '7898')]
    (tax / 'nodes.dmp').write_text(''.join('\t|\t'.join(
        [tid, parent, 'no rank', '', '0', '0', '1', '0', '2', '0', '0', '0', '']) + '\t|\n'
        for tid, parent in nodes))
    (tax / 'names.dmp').write_text(''.join(
        f'{tid}\t|\tSynthetic {tid}\t|\t\t|\tscientific name\t|\n' for tid, _ in nodes))
    table(folder / 'rna.csv', ['readset_id', 'library_id', 'layout', 'fastq_1', 'fastq_2'],
          [['paired', 'L1', 'paired', 'r1.fq', 'r2.fq']])
    header = ['assembly_id', 'assembly_fasta', 'rna_readsets', 'use_protein_evidence']
    table(folder / 'ready.csv', header, [
        ['combined', 'genome.fa', 'paired', 'true'],
        ['rna', 'genome.fa', 'paired', 'false'],
        ['protein', 'genome.fa', '', 'true'],
        ['genome', 'genome.fa', '', 'false']])
    table(folder / 'blocked.csv', header, [
        ['blocked_a', 'genome.fa', 'missing', 'true'],
        ['blocked_b', 'genome.fa', 'missing', 'false'],
        ['independent', 'genome.fa', '', 'true']])
    table(folder / 'bad_rna.csv', ['readset_id', 'library_id', 'layout', 'fastq_1'],
          [['missing', 'L2', 'single', 'absent.fq']])
    # Reproduce a trailing-comma schema error rather than an ordinary empty cell.
    (folder / 'malformed.csv').write_text(','.join(header) + '\nA,genome.fa,,false,\n')
    return tax


def run_case(name, sheet, rna, expected_status, expected_states, folder, tax):
    case = folder / name
    case.mkdir()
    output = case / 'results'
    output.mkdir()
    command = ['nextflow', '-log', str(case / 'nextflow.log'),
        '-C', str(ROOT / 'nextflow.config'), 'run', str(ROOT / 'validation/entrypoints/preflight.nf'),
        '-profile', 'slurm', '-w', str(case / 'work'),
        '--sample_sheet', str(sheet), '--rna_sheet', str(rna),
        '--taxid', '999', '--taxdump_dir', str(tax), '--outdir', str(output),
        '--ref_proteome_accession', 'GCF_024713975.1']
    (case / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
    print(f'Running {name} (tiny synthetic preflight only)...', flush=True)
    with (case / 'console.log').open('w') as log:
        result = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    report = output / 'preflight'
    manifest_path = report / 'resolved_inputs.json'
    if not manifest_path.is_file():
        console_path = case / 'console.log'
        tail = '\n'.join(console_path.read_text(errors='replace').splitlines()[-35:])
        raise RuntimeError(f'{name}: Nextflow exited {result.returncode} without a published manifest. '
                           f'Inspect {console_path} and {case / "nextflow.log"}.\n'
                           f'Last console lines:\n{tail}')
    manifest = json.loads(manifest_path.read_text())
    assert (report / 'README.md').is_file(), f'{name}: published report missing'
    assert manifest['status'] == expected_status, f'{name}: {manifest}'
    assert (result.returncode == 0) == (expected_status == 'ready'), (
        f'{name}: unexpected exit {result.returncode}; inspect {case / "console.log"}')
    states = {row['assembly_id']: row['status'] for row in manifest.get('assemblies', [])}
    assert states == expected_states, f'{name}: unexpected assembly states {states}'
    if name == 'ready':
        assemblies = {row['assembly_id']: row for row in manifest['assemblies']}
        assert {key: row['prediction_mode'] for key, row in assemblies.items()} == {
            'combined': 'rna_protein', 'rna': 'rna', 'protein': 'protein', 'genome': 'ab_initio'}
        assert len(manifest['rna_readsets']) == 1, 'Shared RNA was duplicated'
        assert assemblies['combined']['settings']['ref_proteome_accession'] == ['GCF_024713975.1']
        assert assemblies['rna']['settings']['ref_proteome_accession'] == []
        assert assemblies['combined']['settings']['busco_lineage'] == 'actinopterygii_odb10'
        assert assemblies['combined']['origins']['taxid'] == 'global_parameter'
        assert assemblies['combined']['settings']['genetic_code'] == 1
    elif name == 'blocked':
        assert manifest['rna_readsets'][0]['status'] == 'failed'
        assert 'absent.fq' in manifest['rna_readsets'][0]['reason']
        assert 'missing' in manifest['assemblies'][0]['reason']
    elif name == 'malformed':
        assert 'expected 4 fields, found 5' in manifest['error'], manifest
    print(f'PASS: {name}', flush=True)
    return {'case': name, 'passed': True, 'exit_code': result.returncode,
            'manifest_status': manifest['status'], 'report': str(report)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, help='New, empty shared cluster directory for test artifacts')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O/PYTHONOPTIMIZE; acceptance assertions must be enabled')
    folder = Path(args.outdir).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    if any(folder.iterdir()):
        parser.error('--outdir must be empty; use a new directory to avoid stale results')
    summary = {'status': 'failed', 'checks': []}
    try:
        print('Running Python input-contract tests...', flush=True)
        with (folder / 'unit_tests.log').open('w') as log:
            unit = subprocess.run([sys.executable, '-m', 'unittest', 'discover',
                '-s', str(ROOT / 'tests'), '-p', 'test_annotation_inputs.py', '-v'],
                cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if unit.returncode:
            raise RuntimeError(f'Input tests failed; inspect {folder / "unit_tests.log"}')
        summary['checks'].append({'case': 'unit_tests', 'passed': True})
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        version_text = version.stdout + version.stderr
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version_text):
            raise RuntimeError('Load Nextflow 23.10.1 before running this acceptance suite')
        summary['nextflow_version'] = '23.10.1'
        fixture_dir = folder / 'fixtures'
        tax = fixtures(fixture_dir)
        for name, sheet, rna, status, states in (
            ('ready', 'ready.csv', 'rna.csv', 'ready',
             {x: 'ready' for x in ('combined', 'rna', 'protein', 'genome')}),
            ('blocked', 'blocked.csv', 'bad_rna.csv', 'incomplete',
             {'blocked_a': 'blocked', 'blocked_b': 'blocked', 'independent': 'ready'}),
            ('malformed', 'malformed.csv', 'rna.csv', 'failed', {}),
        ):
            summary['checks'].append(run_case(name, fixture_dir / sheet,
                fixture_dir / rna, status, states, folder, tax))
        summary['status'] = 'passed'
    except (OSError, ValueError, AssertionError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', file=sys.stderr)
    (folder / 'acceptance.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {folder / "acceptance.json"}', flush=True)
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
