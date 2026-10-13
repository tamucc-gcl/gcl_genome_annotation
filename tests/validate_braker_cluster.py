#!/usr/bin/env python3
"""Connected tiny real-tool handoff and BRAKER3 runtime acceptance; no gene training."""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

from validate_preflight_cluster import table
from validate_rna_alignment_cluster import create_inputs

ROOT = Path(__file__).resolve().parents[1]


def fixtures(folder):
    create_inputs(folder)
    protein = folder / 'proteins.fa'
    protein.write_text('>synthetic_protein\n' + 'M' + 'ACDEFGHIKLMNPQRSTVWY' * 10 + '\n')
    genome = ''.join((folder / 'genome.fa').read_text().splitlines()[1:])
    (folder / 'repeats.fa').write_text('>synthetic_repeat#Unknown\n' + genome[:200] + '\n')
    # Explicit local proteins prevent taxonomy defaults from requesting downloads.
    table(folder / 'combined.csv',
          ['assembly_id', 'assembly_fasta', 'rna_readsets', 'use_protein_evidence',
           'orthodb_clade', 'protein_fastas', 'ref_proteome_accession'],
          [['A', 'genome.fa', 'P1;P2;S', 'true', 'none', str(protein), ''],
           ['B', 'genome.fa', 'P1', 'true', 'none', str(protein), '']])
    table(folder / 'reads.csv',
          ['readset_id', 'library_id', 'sample_id', 'layout', 'fastq_1', 'fastq_2',
           'trim_adapters', 'trim_quality'],
          [[rid, lid, lid, layout, str(folder / f'{rid}_R1.fq'),
            str(folder / f'{rid}_R2.fq') if layout == 'paired' else '', 'false', 'false']
           for rid, lid, layout in [('P1', 'P', 'paired'), ('P2', 'P', 'paired'), ('S', 'S', 'single')]])


def run_case(root, inputs, entry, suffix, resume=False):
    case = root / ('connected' if entry == 'main.nf' else 'prediction_only')
    case.mkdir(exist_ok=True)
    output = case / 'results'
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / entry), '-profile', 'slurm', '-w', str(case / 'work'),
               '--sample_sheet', str(inputs / 'combined.csv'), '--rna_sheet', str(inputs / 'reads.csv'),
               '--taxid', '999', '--taxdump_dir', str(inputs / 'taxdump'), '--outdir', str(output),
               '--braker_check_only', 'true', '--braker_cpus', '2', '--braker_threads', '2',
               '--braker_rounds', '1', '--min_contig_length', '0',
               '--repeat_library', str(inputs / 'repeats.fa'), '--n_rm_chunks', '2',
               '--repeat_mask_cpus', '2', '--repeat_mask_memory', '8 GB',
               '--rna_alignment_cpus', '4', '--rna_alignment_threads', '2',
               '--rna_alignment_memory', '8 GB', '--rna_strand_min_informative', '20',
               '--protein_source_store', str(root / 'protein_store'),
               '--orthodb_database', str(root / 'orthodb_store')]
    if entry != 'main.nf':
        command += ['--braker_inputs_root', str(root / 'connected/results')]
    if resume:
        command.append('-resume')
    print(f'Running {entry} ({suffix}; real tools, BRAKER software check only)...', flush=True)
    (case / f'{suffix}.command.json').write_text(json.dumps(command, indent=2))
    with (case / f'{suffix}.console.log').open('w') as log:
        process = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    report = output / 'structural_annotation/prediction_status.json'
    if process.returncode or not report.exists():
        tail = '\n'.join((case / f'{suffix}.console.log').read_text(errors='replace').splitlines()[-60:])
        raise RuntimeError(f'{entry}: exit {process.returncode}, report exists={report.exists()}.\n{tail}')
    data = json.loads(report.read_text())
    assert data['status'] == 'software_checked'
    rows = {r['assembly_id']: r for r in data['assemblies']}
    assert set(rows) == {'A', 'B'}
    assert [r['readset_id'] for r in rows['A']['readsets']] == ['P1', 'P2', 'S']
    assert [r['readset_id'] for r in rows['B']['readsets']] == ['P1']
    for aid, count in [('A', 3), ('B', 1)]:
        runtime = output / 'structural_annotation/braker_runtime' / aid
        provenance = json.loads((runtime / 'inputs.json').read_text())
        assert provenance['execution_mode'] == 'software_check'
        assert re.search(r'\bversion\s+3\.', (runtime / 'braker_version.txt').read_text())
        args = json.loads((runtime / 'command.json').read_text())
        bam_option = next(a for a in args if a.startswith('--bam='))
        assert len(bam_option.split('=', 1)[1].split(',')) == count
        assert '--checkSoftware' in args and not any(a.startswith('--stranded') for a in args)
        handoff = json.loads((output / 'structural_annotation/braker_inputs' / f'{aid}.inputs.json').read_text())
        assert handoff['status'] == 'ready'
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    # Preserve each run's trace before the next run overwrites the shared report.
    (case / f'{suffix}.trace.json').write_text(json.dumps(trace, indent=2) + '\n')
    cache_details = '; '.join(f"{r['name']}: {r['status']} ({r.get('hash', '?')})" for r in trace
                             if r['name'].startswith(('FILTER_ASSEMBLY', 'COMBINE_PROTEIN_EVIDENCE',
                                                      'INDEX_RNA_REFERENCE', 'ALIGN_RNA_READSET',
                                                      'VALIDATE_BRAKER_INPUTS', 'CHECK_BRAKER_SOFTWARE')))
    assert not any(r['name'].split(' (')[0] == 'BRAKER' for r in trace), 'Unexpected training task'
    for name in ('VALIDATE_BRAKER_INPUTS', 'CHECK_BRAKER_SOFTWARE'):
        tasks = [r for r in trace if r['name'].startswith(name)]
        assert len(tasks) == 2
        if resume and suffix != 'recover':
            assert all(r['status'] == 'CACHED' for r in tasks), f'{name} not cached. Task cache details: {cache_details}'
    if entry == 'main.nf':
        alignments = [r for r in trace if r['name'].startswith('ALIGN_RNA_READSET')]
        assert len(alignments) == 4
        if resume and suffix != 'recover':
            assert all(r['status'] == 'CACHED' for r in alignments), f'Alignment unexpectedly reran. Task cache details: {cache_details}'
    print(f'PASS: {entry} ({suffix})', flush=True)
    return dict(entry=entry, run=suffix, passed=True, exit_code=process.returncode)


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
    inputs = root / 'fixtures'
    if args.resume and not (inputs / 'combined.csv').exists():
        parser.error('--resume requires the original fixtures')
    summary = dict(status='failed', checks=[], scope='Connected preparation/alignment/handoff and actual BRAKER3 software checks; no training or prediction')
    try:
        with (root / 'unit_tests.log').open('w') as log:
            unit = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'), '-p', 'test_*.py'],
                                  cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if unit.returncode:
            raise RuntimeError('Python tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests', passed=True))
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version.stdout + version.stderr):
            raise RuntimeError('Load Nextflow 23.10.1')
        if not args.resume:
            fixtures(inputs)
        summary['checks'].append(run_case(root, inputs, 'main.nf', 'recover' if args.resume else 'fresh', args.resume))
        summary['checks'].append(run_case(root, inputs, 'main.nf', 'resume', True))
        summary['checks'].append(run_case(root, inputs, 'braker_prediction.nf', 'fresh', args.resume))
        summary['checks'].append(run_case(root, inputs, 'braker_prediction.nf', 'resume', True))
        summary['status'] = 'passed'
    except (AssertionError, ValueError, OSError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', flush=True)
    (root / 'acceptance.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {root / "acceptance.json"}')
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
