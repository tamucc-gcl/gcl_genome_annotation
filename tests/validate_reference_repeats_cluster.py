#!/usr/bin/env python3
"""Cluster acceptance with actual RepeatMasker: integrity, unsplit comparison and resume.

The supplied synthetic library deliberately bypasses discovery/classification/database
filtering. Those branches require the subsequent real-assembly pilot.
"""
import argparse
import csv
import json
import random
import re
import subprocess
import sys
from pathlib import Path
from validate_preflight_cluster import fixtures, table

ROOT = Path(__file__).resolve().parents[1]


def make_inputs(folder):
    fixtures(folder)
    rng = random.Random(42)
    sequence = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = sequence(600)
    (folder / 'library.fa').write_text('>fixture%name#DNA/hAT\n' + repeat + '\n')
    # Whole contigs, long IDs, internal double underscores and an apparent range suffix.
    (folder / 'genome.fa').write_text('>long__' + 'x' * 75 + ' original description\n' +
                                     sequence(400) + repeat + sequence(400) + '\n' +
                                     '>contig_100-500\n' + sequence(300) + repeat + sequence(300) + '\n' +
                                     '>small\n' + sequence(200) + '\n')
    (folder / 'invalid.fa').write_text('>collision|id\nACGT\n>collision_id\n' + sequence(1200) + '\n')
    (folder / 'invalid_library.fa').write_text('>duplicate\nACGT\n>duplicate\nTGCA\n')
    header = ['assembly_id', 'assembly_fasta', 'use_protein_evidence', 'min_contig_length']
    ready = [['A', 'genome.fa', 'false', '1000'], ['B', 'genome.fa', 'false', '0']]
    table(folder / 'ready.csv', header, ready)
    table(folder / 'added.csv', header, ready + [['C', 'genome.fa', 'false', '1000']])
    table(folder / 'blocked.csv', header, [ready[0], ['bad', 'invalid.fa', 'false', '1000']])


def run_case(root, inputs, name, sheet, expected, suffix='fresh', resume=False, library='library.fa'):
    case = root / name
    case.mkdir(exist_ok=True)
    output = case / 'results'
    output.mkdir(exist_ok=True)
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / 'validation/entrypoints/reference_repeats.nf'), '-profile', 'slurm', '-w', str(case / 'work'),
               '--sample_sheet', str(inputs / sheet), '--taxid', '999',
               '--taxdump_dir', str(inputs / 'taxdump'), '--outdir', str(output),
               '--repeat_library', str(inputs / library), '--repeat_mask_benchmark', 'true',
               '--n_rm_chunks', '8', '--repeat_mask_cpus', '4',
               '--repeat_mask_memory', '8 GB', '--repeat_mask_time', '2h', '--publish_dir_mode', 'copy']
    if resume:
        command.append('-resume')
    (case / f'{suffix}.command.json').write_text(json.dumps(command, indent=2) + '\n')
    print(f'Running reference/repeats {name} ({suffix}; actual RepeatMasker)...', flush=True)
    console = case / f'{suffix}.console.log'
    with console.open('w') as log:
        run = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    report = output / 'reference_repeats/reference_repeat_status.json'
    if not report.is_file():
        tail = '\n'.join(console.read_text(errors='replace').splitlines()[-40:])
        raise RuntimeError(f'{name}: exit {run.returncode}; missing reference/repeat report.\n{tail}')
    data = json.loads(report.read_text())
    assert {row['assembly_id']: row['status'] for row in data['assemblies']} == expected, data
    assert (run.returncode == 0) == all(value == 'ready' for value in expected.values()), console
    for row in data['assemblies']:
        aid = row['assembly_id']
        if row['status'] != 'ready':
            continue
        mask = row['repeats']
        assert mask['repeat_interval_union_bases'] > 0, 'Fixture failed to exercise real masking'
        assert mask['lowercase_acgt_bases'] > 0
        assert mask['contigs'] == (3 if aid == 'B' else 2)
        accounting = json.loads((output / f'repeats/batches/{aid}.batches.json').read_text())
        assert accounting['actual_batches'] == mask['contigs']
        assert accounting['requested_batches'] == 8
        assert sum(accounting['batch_bases']) == mask['reference_bases']
        benchmark = json.loads((output / f'reference_repeats/benchmarks/{aid}.benchmark.json').read_text())
        assert benchmark['status'] == 'ready', benchmark
        gff = (output / f'repeats/{aid}.repeat.gff3').read_text()
        ids = re.findall(r'\bID=([^;\s]+)', gff)
        assert ids and len(ids) == len(set(ids)), 'GFF feature IDs must be unique'
        assert 'repeat_linkage_id=' in gff
        assert 'fixture%25name' in gff, 'GFF attributes did not escape the repeat identifier'
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    forbidden = ('REPEAT_MODEL', 'CLASSIFY_REPEATS', 'FILTER_REPEAT_LIBRARY', 'CURATE_DFAM')
    assert not any(any(task in row['name'] for task in forbidden) for row in trace), 'Supplied-library route ran discovery'
    tasks = ('FILTER_ASSEMBLY', 'PREP_RM_BATCHES', 'REPEAT_MASK_BATCH', 'MERGE_RM_RESULTS', 'VALIDATE_REFERENCE_REPEATS')
    if resume:
        for task in tasks:
            rows = [row for row in trace if re.search(r'(?:^|:)' + task + r' \((?:A|B)(?:__benchmark)?(?::|\))', row['name'])]
            assert rows, f'Missing trace rows for {task}'
            assert all(row['status'] == 'CACHED' for row in rows), f'{task}: unchanged assemblies reran'
        if suffix == 'added':
            rows = [row for row in trace if 'REPEAT_MASK_BATCH (C:' in row['name']]
            assert rows and all(row['status'] == 'COMPLETED' for row in rows), 'Added assembly did not run'
    if name == 'blocked':
        assert not any('REPEAT_MASK_BATCH (bad:' in row['name'] for row in trace), 'Invalid reference was masked'
    if name == 'invalid_library':
        assert not any('REPEAT_MASK_BATCH' in row['name'] for row in trace), 'Invalid library reached masking'
    print(f'PASS: reference/repeats {name} ({suffix})', flush=True)
    return {'case': name, 'run': suffix, 'passed': True, 'exit_code': run.returncode}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O/PYTHONOPTIMIZE')
    root = Path(args.outdir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        parser.error('--outdir must be empty')
    summary = {'status': 'failed', 'checks': [], 'scope': 'reference integrity and actual custom-library masking; not de novo discovery'}
    try:
        print('Running Python integrity-contract tests...', flush=True)
        with (root / 'unit_tests.log').open('w') as log:
            tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'),
                                    '-p', 'test_*.py', '-v'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        assert tests.returncode == 0, f'Unit tests failed; inspect {root / "unit_tests.log"}'
        summary['checks'].append({'case': 'unit_tests', 'passed': True})
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True, check=True)
        match = re.search(r'version\s+([\d.]+)', version.stdout + version.stderr)
        summary['nextflow_version'] = match.group(1) if match else 'unknown'
        inputs = root / 'inputs'
        make_inputs(inputs)
        ready = {'A': 'ready', 'B': 'ready'}
        for suffix, sheet, states, resume in (('fresh', 'ready.csv', ready, False),
                                             ('resume', 'ready.csv', ready, True),
                                             ('added', 'added.csv', dict(ready, C='ready'), True)):
            summary['checks'].append(run_case(root, inputs, 'ready', sheet, states, suffix, resume))
        summary['checks'].append(run_case(root, inputs, 'blocked', 'blocked.csv', {'A': 'ready', 'bad': 'blocked'}))
        summary['checks'].append(run_case(root, inputs, 'invalid_library', 'ready.csv',
                                         {'A': 'blocked', 'B': 'blocked'}, library='invalid_library.fa'))
        summary['status'] = 'passed'
    except (AssertionError, OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        summary['error'] = str(exc)
        print('FAIL: ' + str(exc), flush=True)
    path = root / 'acceptance.json'
    path.write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {path}', flush=True)
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
