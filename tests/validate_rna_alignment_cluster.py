#!/usr/bin/env python3
"""Actual HISAT2/BAM acceptance: splices, mates, lane grouping, resume and strand-only changes."""
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
sys.path.insert(0, str(ROOT / 'py_scripts'))
from prepare_rna import inspect
from reference_repeats import reference


def create_inputs(root):
    fixtures(root)
    rng = random.Random(9473)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    exons = [dna(350) for _ in range(3)]
    genome = dna(500) + exons[0] + 'GT' + dna(116) + 'AG' + exons[1] + 'GT' + dna(116) + 'AG' + exons[2] + dna(500)
    transcript = ''.join(exons)
    (root / 'genome.fa').write_text('>chr1\n' + genome + '\n')
    (root / 'different.fa').write_text('>chr2\n' + genome + '\n')
    rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    selected = []
    for rid, lid, layout in (('P1', 'P', 'paired'), ('P2', 'P', 'paired'), ('S', 'S', 'single')):
        paths = [root / f'{rid}_R{i}.fq' for i in range(1, 3 if layout == 'paired' else 2)]
        records = [[] for _ in paths]
        count = 0
        for offset in range(0, 800, 5):
            for repeat in range(4):
                count += 1
                seqs = [transcript[offset:offset + 100], rc(transcript[offset + 150:offset + 250])] if layout == 'paired' else [rc(transcript[offset:offset + 100])]
                for mate, seq in enumerate(seqs, 1):
                    records[mate - 1].append(f'@{rid}_{count}/{mate}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
        for path, rows in zip(paths, records):
            path.write_text(''.join(rows))
        selected.append(dict(readset_id=rid, library_id=lid, sample_id=lid, layout=layout,
                             status='ready', strandedness='auto', read_paths=list(map(str, paths)),
                             **inspect(paths, layout)))
    prepared = root / 'prepared.json'
    prepared.write_text(json.dumps(dict(status='ready', readsets=selected)))
    table(root / 'reads.csv', ['readset_id', 'library_id', 'sample_id', 'layout', 'fastq_1', 'fastq_2'],
          [[row['readset_id'], row['library_id'], row['sample_id'], row['layout'], row['read_paths'][0],
            row['read_paths'][1] if row['layout'] == 'paired' else ''] for row in selected])
    refs = root / 'references'
    for aid in ('A', 'B'):
        result = reference(root / 'genome.fa', aid, 0, refs / aid)
        assert result['status'] == 'ready'
    assert reference(root / 'different.fa', 'C', 0, refs / 'C')['status'] == 'ready'
    header = ['assembly_id', 'assembly_fasta', 'rna_readsets', 'use_protein_evidence']
    rows = [['A', 'genome.fa', 'P1;P2;S', 'false'], ['B', 'genome.fa', 'P1', 'false']]
    table(root / 'ready.csv', header, rows)
    table(root / 'added.csv', header, rows + [['C', 'different.fa', 'S', 'false']])
    table(root / 'none.csv', header, [['N', 'genome.fa', '', 'false']])


def run_case(root, inputs, name, sheet, expected, suffix='fresh', resume=False, strand=None, minimum=20):
    case = root / name
    case.mkdir(exist_ok=True)
    output = case / 'results'
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / 'rna_alignment.nf'), '-profile', 'slurm', '-w', str(case / 'work'),
               '--sample_sheet', str(inputs / sheet), '--rna_sheet', str(inputs / 'reads.csv'),
               '--rna_prepared_status', str(inputs / 'prepared.json'), '--rna_reference_dir', str(inputs / 'references'),
               '--taxid', '999', '--taxdump_dir', str(inputs / 'taxdump'), '--outdir', str(output),
               '--rna_alignment_memory', '8 GB',
               '--rna_alignment_cpus', '4', '--rna_alignment_threads', '2',
               '--rna_strand_min_informative', str(minimum)]
    if resume:
        command.append('-resume')
    if strand:
        command += ['--strandedness', strand]
    print(f'Running RNA alignment {name} ({suffix}; actual HISAT2)...', flush=True)
    (case / f'{suffix}.command.json').write_text(json.dumps(command, indent=2))
    with (case / f'{suffix}.console.log').open('w') as log:
        process = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    status = output / 'rna_alignment/rna_alignment_status.json'
    if not status.exists():
        tail = '\n'.join((case / f'{suffix}.console.log').read_text(errors='replace').splitlines()[-40:])
        raise RuntimeError(f'{name}: exit {process.returncode}; missing alignment report.\n{tail}')
    data = json.loads(status.read_text())
    alignment_dir = output / 'rna_alignment'
    assert '(multiqc/multiqc_report.html)' in (alignment_dir / 'README.md').read_text()
    assert (alignment_dir / 'multiqc/multiqc_report.html').is_file()
    if data['readsets']:
        expected_samples = {f"{row['assembly_id']}/{row['readset_id']}" for row in data['readsets']}
        for module in ('bowtie2', 'samtools_stats', 'samtools_flagstat'):
            with (alignment_dir / f'multiqc/multiqc_report_data/multiqc_{module}.txt').open() as handle:
                module_rows = list(csv.DictReader(handle, delimiter='\t'))
            assert {next(iter(row.values())) for row in module_rows} == expected_samples, f'MultiQC {module}: lost/mislabelled sample'
    else:
        assert 'No RNA alignment was requested' in (alignment_dir / 'multiqc/multiqc_report.html').read_text()
    assert {row['assembly_id']: row['status'] for row in data['assemblies']} == expected
    assert (process.returncode == 0) == all(value != 'blocked' for value in expected.values())
    with (output / 'pipeline/pipeline_trace.txt').open() as handle:
        trace = list(csv.DictReader(handle, delimiter='\t'))
    if name == 'ready':
        expected_counts = {'INDEX_RNA_REFERENCE': 2 if suffix == 'added' else 1,
                           'ALIGN_RNA_READSET': 5 if suffix == 'added' else 4}
        for task, count in expected_counts.items():
            jobs = [row for row in trace if row['name'].startswith(task)]
            assert len(jobs) == count, f'{task}: missing/duplicate tasks'
            if resume and suffix != 'recover':
                cached = count - (1 if suffix == 'added' and task in ('INDEX_RNA_REFERENCE', 'ALIGN_RNA_READSET') else 0)
                assert sum(row['status'] == 'CACHED' for row in jobs) == cached, f'{task}: unexpected invalidation'
        if suffix == 'strand_change':
            strand_jobs = [row for row in trace if row['name'].startswith('RESOLVE_RNA_STRAND')]
            assert all(row['status'] == 'COMPLETED' for row in strand_jobs), 'Strand change not applied'
        for row in data['readsets']:
            assert row['mapped_fraction'] > .95
            assert row['junctions'] == 2, 'Known synthetic splice junctions missing/extra'
            assert row['strand']['strandedness'] == (strand or ('reverse' if row['readset_id'] == 'S' else 'forward'))
        libs = {(row['assembly_id'], row['library_id']): row for row in data['libraries']}
        assert libs[('A', 'P')]['readsets'] == ['P1', 'P2']
        assert libs[('A', 'S')]['readsets'] == ['S']
        assert libs[('B', 'P')]['readsets'] == ['P1']
        assert len(libs[('A', 'P')]['bams']) == 2, 'Technical lanes were merged/dropped'
        for row in data['libraries']:
            for path in row['bams'] + row['bam_indices']:
                assert (output / 'rna_alignment' / path).is_file()
        assert not any(row['name'].startswith('MERGE_RNA') for row in trace)
    if name == 'blocked':
        assert all(row['strand']['status'] == 'blocked' for row in data['readsets'])
    if name == 'none':
        assert not any(row['name'].startswith(('INDEX_RNA_REFERENCE', 'ALIGN_RNA_READSET', 'MERGE_RNA_LIBRARY')) for row in trace)
    print(f'PASS: RNA alignment {name} ({suffix})', flush=True)
    return dict(case=name, run=suffix, passed=True, exit_code=process.returncode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--resume', action='store_true', help='Reuse existing fixtures and resume an interrupted ready run')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O/PYTHONOPTIMIZE')
    root = Path(args.outdir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()) and not args.resume:
        parser.error('--outdir must be empty')
    if args.resume and not all((root / path).exists() for path in
                               ('fixtures/prepared.json', 'fixtures/ready.csv', 'ready/.nextflow')):
        parser.error('--resume needs the original fixtures and ready/.nextflow cache')
    summary = dict(status='failed', checks=[], scope='Actual HISAT2 alignment and BAM handoff; no repeats or gene prediction')
    try:
        with (root / 'unit_tests.log').open('w') as log:
            unit = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'), '-p', 'test_*.py'],
                                  cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        if unit.returncode:
            raise RuntimeError('Python contract tests failed; inspect unit_tests.log')
        summary['checks'].append(dict(case='unit_tests', passed=True))
        version = subprocess.run(['nextflow', '-version'], capture_output=True, text=True)
        if version.returncode or not re.search(r'\bversion\s+23\.10\.1\b', version.stdout + version.stderr):
            raise RuntimeError('Load Nextflow 23.10.1')
        inputs = root / 'fixtures'
        if not args.resume:
            create_inputs(inputs)
        expected = dict(A='ready', B='ready')
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready.csv', expected,
                                         'recover' if args.resume else 'fresh', resume=args.resume))
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready.csv', expected, 'resume', resume=True))
        summary['checks'].append(run_case(root, inputs, 'ready', 'ready.csv', expected, 'strand_change', resume=True, strand='reverse'))
        summary['checks'].append(run_case(root, inputs, 'ready', 'added.csv', {**expected, 'C': 'ready'}, 'added', resume=True))
        summary['checks'].append(run_case(root, inputs, 'blocked', 'ready.csv', dict(A='blocked', B='blocked'), minimum=1000000))
        summary['checks'].append(run_case(root, inputs, 'none', 'none.csv', dict(N='not_requested')))
        summary['status'] = 'passed'
    except (AssertionError, ValueError, OSError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', flush=True)
    (root / 'acceptance.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Acceptance summary: {root / "acceptance.json"}')
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
