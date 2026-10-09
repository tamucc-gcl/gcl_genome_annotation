#!/usr/bin/env python3
"""Cluster acceptance for paired/single short reads, failure isolation, QC and resume."""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from validate_preflight_cluster import fixtures, table

ROOT = Path(__file__).resolve().parents[1]


def run_case(root, inputs, name, sheet, expected, resume=False, fraction=None, passthrough=False):
    case = root / name
    case.mkdir(exist_ok=True)
    output = case / 'results'
    output.mkdir(exist_ok=True)
    suffix = 'subset_change' if fraction else ('resume' if resume else 'fresh')
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / 'rna_preparation.nf'), '-profile', 'slurm', '-w', str(case / 'work'),
               '--sample_sheet', str(inputs / sheet), '--rna_sheet', str(inputs / 'reads.csv'),
               '--taxid', '999', '--taxdump_dir', str(inputs / 'taxdump'), '--outdir', str(output)]
    if resume:
        command.append('-resume')
    if fraction:
        command += ['--subset_fraction', str(fraction)]
    if passthrough:
        command += ['--trim_adapters', 'false', '--trim_quality', 'false']
    print(f'Running RNA {name} ({suffix})...', flush=True)
    with (case / f'{suffix}.console.log').open('w') as log:
        process = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    path = output / 'rna/rna_status.json'
    if not path.exists():
        tail = '\n'.join((case / f'{suffix}.console.log').read_text(errors='replace').splitlines()[-40:])
        raise RuntimeError(f'{name}: exit {process.returncode}; missing RNA report.\n{tail}')
    data = json.loads(path.read_text())
    assert {row['assembly_id']: row['status'] for row in data['assemblies']} == expected
    assert (process.returncode == 0) == all(value != 'blocked' for value in expected.values())
    for report_stage, read_stage in (('raw', 'raw'), ('filtered', 'processed')):
        report_dir = output / f'rna/multiqc/{report_stage}'
        assert (report_dir / 'multiqc_report.html').is_file()
        ready = [row for row in data['readsets'] if row['status'] == 'ready']
        if ready:
            with (report_dir / 'multiqc_data/multiqc_fastqc.txt').open() as handle:
                samples = {row['Sample'] for row in csv.DictReader(handle, delimiter='\t')}
            for row in ready:
                expected_samples = {f"{row['readset_id']}_{read_stage}_R{mate}" for mate in
                                    range(1, 3 if row['layout'] == 'paired' else 2)}
                assert expected_samples <= samples, f'{report_stage}: missing FastQC samples'
            other_stage = 'processed' if read_stage == 'raw' else 'raw'
            assert not any(f'_{other_stage}_R' in name for name in samples), 'QC stages mixed'
            if report_stage == 'filtered' and not passthrough:
                assert (report_dir / 'multiqc_data/multiqc_fastp.txt').is_file()
    for row in data['readsets']:
        if row['status'] == 'ready':
            assert row['read_units'] > 0
            assert len(row['read_paths']) == (2 if row['layout'] == 'paired' else 1)
            assert all((output / 'rna' / name).is_file() for name in row['read_paths'])
            rid = row['readset_id']
            # Local validation and all-read selection must publish metadata only.
            assert not list((output / f'rna/raw/{rid}').glob('*.f*q*'))
            if not fraction:
                assert not list((output / f'rna/selected/{rid}').glob('*.f*q*'))
            if passthrough:
                assert not list((output / f'rna/processed/{rid}').glob('*.f*q*'))
                request = json.loads((output / f'rna/raw/{rid}/rna.json').read_text())['request']
                assert row['read_paths'] == [request[key] for key in ('fastq_1', 'fastq_2') if request.get(key)]
            for stage in ('raw', 'processed'):
                assert list((output / f'rna/qc/{stage}/{rid}_{stage}').glob('*_fastqc.zip'))
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    if name in ('ready', 'passthrough'):
        for task in ('ACQUIRE_RNA_READSET', 'QC_RNA_RAW', 'TRIM_RNA_READSET', 'QC_RNA_PROCESSED', 'SUBSET_RNA_READSET'):
            jobs = [row for row in trace if row['name'].startswith(task)]
            assert len(jobs) == 2, f'{task}: shared readsets were duplicated/missing'
            if resume and not (fraction and task == 'SUBSET_RNA_READSET'):
                assert all(row['status'] == 'CACHED' for row in jobs), f'{task}: unexpected rerun'
            if fraction and task == 'SUBSET_RNA_READSET':
                assert all(row['status'] == 'COMPLETED' for row in jobs), 'Subset settings did not invalidate selection'
    print(f'PASS: RNA {name} ({suffix})', flush=True)
    return {'case': name, 'resume': resume, 'subset_fraction': fraction, 'passed': True, 'exit_code': process.returncode}


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
    summary = {'status': 'failed', 'checks': []}
    try:
        with (root / 'unit_tests.log').open('w') as log:
            tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'),
                                    '-p', 'test_*.py', '-v'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if tests.returncode:
            raise RuntimeError(f'Unit tests failed; inspect {root / "unit_tests.log"}')
        summary['checks'].append({'case': 'unit_tests', 'passed': True})
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version.stdout + version.stderr):
            raise RuntimeError('Load Nextflow 23.10.1')
        inputs = root / 'fixtures'
        fixtures(inputs)
        sequence = 'ACGT' * 25
        for name, mate in (('pair1', 1), ('pair2', 2), ('single', 1)):
            (inputs / f'{name}.fq').write_text(''.join(
                f'@read{index}/{mate}\n{sequence}\n+\n' + 'I' * 100 + '\n' for index in range(256)))
        (inputs / 'bad2.fq').write_text('@wrong/2\nACGT\n+\nIIII\n')
        table(inputs / 'reads.csv', ['readset_id', 'library_id', 'layout', 'fastq_1', 'fastq_2', 'strandedness'], [
            ['P', 'paired_lib', 'paired', 'pair1.fq', 'pair2.fq', 'auto'],
            ['S', 'single_lib', 'single', 'single.fq', '', 'reverse'],
            ['bad', 'bad_lib', 'paired', 'pair1.fq', 'bad2.fq', 'auto']])
        # A source cannot be assigned twice: create a separate copy for the bad row.
        (inputs / 'bad1.fq').write_text((inputs / 'pair1.fq').read_text())
        text = (inputs / 'reads.csv').read_text().replace('bad_lib,paired,pair1.fq', 'bad_lib,paired,bad1.fq')
        (inputs / 'reads.csv').write_text(text)
        header = ['assembly_id', 'assembly_fasta', 'rna_readsets', 'use_protein_evidence']
        table(inputs / 'ready_rna.csv', header, [
            ['A', 'genome.fa', 'P;S', 'false'], ['B', 'genome.fa', 'P', 'false'], ['none', 'genome.fa', '', 'false']])
        table(inputs / 'blocked_rna.csv', header, [
            ['A', 'genome.fa', 'bad', 'false'], ['B', 'genome.fa', 'bad', 'false'], ['C', 'genome.fa', 'S', 'false']])
        table(inputs / 'no_rna.csv', header, [['none', 'genome.fa', '', 'false']])
        expected = {'A': 'ready', 'B': 'ready', 'none': 'not_requested'}
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready_rna.csv', expected))
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready_rna.csv', expected, resume=True))
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready_rna.csv', expected, resume=True, fraction=0.5))
        summary['checks'].append(run_case(root, inputs, 'passthrough', 'ready_rna.csv', expected, passthrough=True))
        summary['checks'].append(run_case(root, inputs, 'passthrough', 'ready_rna.csv', expected, resume=True, passthrough=True))
        summary['checks'].append(run_case(root, inputs, 'blocked', 'blocked_rna.csv',
                                         {'A': 'blocked', 'B': 'blocked', 'C': 'ready'}))
        summary['checks'].append(run_case(root, inputs, 'none', 'no_rna.csv', {'none': 'not_requested'}))
        summary['status'] = 'passed'
    except (OSError, ValueError, RuntimeError, AssertionError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', file=sys.stderr)
    path = root / 'acceptance.json'
    path.write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {path}')
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
