#!/usr/bin/env python3
"""Tiny real-miniprot AED pilot: four modes, retained GTF reuse and resume.

No BRAKER, mapping, repeats or downloads. Standalone StringTie fallback and
full-database performance require separate acceptance on real BAMs/data.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'py_scripts'))
from prepare_protein_evidence import checksum


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def fixtures(root):
    accepted = root / 'accepted'
    # Diverse coding sequence, avoiding a low-complexity all-methionine fixture.
    codons = dict(zip('ACDEFGHIKLMNPQRSTVWY',
                      ['GCT','TGT','GAT','GAA','TTT','GGT','CAT','ATT','AAA','CTG',
                       'ATG','AAT','CCT','CAA','CGT','TCT','ACT','GTT','TGG','TAT']))
    rng = random.Random(42)
    protein = 'M' + ''.join(rng.choice(list(codons)) for _ in range(299))
    coding = ''.join(codons[a] for a in protein) + 'TAA'
    genome = ''.join(rng.choice('ACGT') for _ in range(500)) + coding + ''.join(rng.choice('ACGT') for _ in range(2000))
    contigs = [dict(id='chr1', length=len(genome), uppercase_sha256=hashlib.sha256(genome.encode()).hexdigest())]
    identity = dict(contigs=contigs, reference_identity=hashlib.sha256(
        json.dumps(contigs, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
    predictions, assemblies = [], []
    for aid, mode in [('A', 'rna_protein'), ('B', 'rna'), ('C', 'protein'), ('D', 'ab_initio')]:
        assemblies.append(dict(assembly_id=aid, prediction_mode=mode))
        gff = accepted / f'structural_annotation/braker/{aid}/{aid}.gff3'
        gff.parent.mkdir(parents=True)
        gff.write_text('##gff-version 3\n' + '\n'.join(
            f'chr1\tBRAKER\t{kind}\t501\t1403\t.\t+\t{phase}\t{attrs}'
            for kind, phase, attrs in [('gene','.','ID=g1'), ('mRNA','.','ID=t1;Parent=g1'),
                                        ('exon','.','Parent=t1'), ('CDS','0','Parent=t1')]) + '\n')
        row = dict(assembly_id=aid, mode=mode, status='ready', reference_identity=identity['reference_identity'],
                   gff3_sha256=checksum(gff))
        predictions.append(row)
        write_json(accepted / f'structural_annotation/braker_status/{aid}.prediction.json', row)
        write_json(accepted / f'filtered_assemblies/{aid}/{aid}.contigs.json', identity)
        genome_path = accepted / f'repeats/{aid}.masked.fasta'
        genome_path.parent.mkdir(exist_ok=True)
        genome_path.write_text('>chr1\n' + genome + '\n')
        if mode in ('rna_protein', 'protein'):
            p = accepted / f'protein_evidence/assemblies/{aid}/proteins.fa'
            p.parent.mkdir(parents=True)
            p.write_text('>reference_protein\n' + protein + '\n')
            write_json(p.parent / 'evidence.json', dict(status='ready', proteins_sha256=checksum(p)))
        if mode in ('rna_protein', 'rna'):
            r = gff.parent / 'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
            r.parent.mkdir(parents=True)
            r.write_text('# StringTie-shaped retained GTF fixture, not a measured transcript assembly\n' +
                         '\n'.join(f'chr1\tStringTie\t{kind}\t501\t1403\t.\t+\t.\tgene_id "rna1"; transcript_id "rna1.t1";'
                                   for kind in ('transcript', 'exon')) + '\n')
    write_json(accepted / 'preflight/resolved_inputs.json', dict(status='ready', assemblies=assemblies))
    write_json(accepted / 'structural_annotation/prediction_status.json', dict(status='ready', assemblies=predictions))


def run_case(root, suffix, resume=False, blocked=False):
    case = root / 'pilot'
    case.mkdir(exist_ok=True)
    results = case / 'results'
    command = ['nextflow', '-log', str(case / f'{suffix}.nextflow.log'), '-C', str(ROOT / 'nextflow.config'),
               'run', str(ROOT / 'evidence_scoring.nf'), '-profile', 'slurm', '-w', str(case / 'work'),
               '--aed_inputs_root', str(root / 'accepted'), '--outdir', str(results),
               '--aed_alignment_cpus', '2', '--aed_alignment_memory', '8 GB', '--aed_alignment_time', '1h']
    if resume:
        command.append('-resume')
    write_json(case / f'{suffix}.command.json', command)
    print(f'Running evidence AED ({suffix}; actual miniprot)...', flush=True)
    with (case / f'{suffix}.console.log').open('w') as log:
        result = subprocess.run(command, cwd=case, stdout=log, stderr=subprocess.STDOUT)
    console = (case / f'{suffix}.console.log').read_text(errors='replace')
    if blocked:
        if not result.returncode or 'Models/reference differ from accepted BRAKER validation' not in console:
            raise RuntimeError('Changed model was not blocked as expected.\n' + console[-8000:])
    else:
        if result.returncode:
            raise RuntimeError(f'{suffix}: exit {result.returncode}\n' + console[-8000:])
        for aid in 'ABCD':
            folder = results / f'structural_annotation/evidence_aed/{aid}'
            status = json.loads((folder / 'aed_status.json').read_text())
            assert status['status'] == 'ready' and status['genes'] == status['transcripts'] == 1
            with (folder / 'aed_scores.tsv').open() as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            assert status['evidence']['rna']['available'] == (aid in 'AB')
            assert status['evidence']['protein']['available'] == (aid in 'AC')
            if aid in 'AB':
                assert float(row['aed']) == 0
            elif aid == 'C':
                assert 0 <= float(row['aed']) < .05 and status['evidence']['protein']['hits'] > 0
            else:
                assert row['aed'] == '.' and row['status'] == 'unavailable'
            assert (folder / 'README.md').is_file()
        with (results / 'pipeline/pipeline_trace.txt').open() as handle:
            tasks = list(csv.DictReader(handle, delimiter='\t'))
        assert not any(re.search(r'(BRAKER|ASSEMBLE_AED_RNA|MERGE_AED_RNA|REPEAT_MODEL|ALIGN_RNA_READSET)', t['name']) for t in tasks)
        if suffix == 'resume':
            expected = {'ALIGN_AED_PROTEINS': 2, 'NORMALIZE_AED_EVIDENCE': 4, 'SCORE_EVIDENCE_AED': 4}
            for name, count in expected.items():
                found = [t for t in tasks if t['name'].startswith(name + ' (')]
                assert len(found) == count and all(t['status'] == 'CACHED' for t in found), f'{name} not cached: {found}'
    print(f'PASS: evidence AED {suffix}', flush=True)
    return dict(case=suffix, passed=True, exit_code=result.returncode)


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
    if args.resume and not (root / 'accepted/preflight/resolved_inputs.json').is_file():
        parser.error('--resume requires the original fixtures')
    summary = dict(status='failed', checks=[], scope='Real miniprot, retained GTF, four evidence modes, cache and changed-model blocking; standalone StringTie and full database pilot pending')
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
        summary['checks'].append(run_case(root, 'recover' if args.resume else 'fresh', args.resume))
        summary['checks'].append(run_case(root, 'resume', True))
        gff = root / 'accepted/structural_annotation/braker/A/A.gff3'
        original = gff.read_bytes()
        try:
            gff.write_bytes(original + b'# changed after validation\n')
            summary['checks'].append(run_case(root, 'changed_model', True, True))
        finally:
            gff.write_bytes(original)
        summary['status'] = 'passed'
    except (AssertionError, ValueError, OSError, RuntimeError) as exc:
        summary['error'] = str(exc)
        print(f'FAIL: {exc}', flush=True)
    write_json(root / 'acceptance.json', summary)
    print(f'Acceptance summary: {root / "acceptance.json"}')
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
