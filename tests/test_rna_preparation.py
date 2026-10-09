"""Cluster unit tests for FASTQ integrity and RNA processing contracts."""
import gzip
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('rna', ROOT / 'py_scripts/prepare_rna.py')
rna = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rna)
spec = importlib.util.spec_from_file_location('rna_report', ROOT / 'py_scripts/report_rna_preparation.py')
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


class RnaTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.inputs = []
        for mate in (1, 2):
            path = self.base / f'input{mate}.fq.gz'
            with gzip.open(path, 'wt') as handle:
                for index in range(100):
                    handle.write(f'@read{index}/{mate}\nACGT\n+\nIIII\n')
            self.inputs.append(path)
        self.request = {'readset_id': 'R', 'layout': 'paired', 'run_accession': None}

    def raw(self, name='raw'):
        path = self.base / name
        result = rna.acquire(self.request, self.inputs, path, 1)
        self.assertEqual(result['status'], 'ready')
        return path

    def test_paired_gzip_and_single_reads(self):
        paired = rna.inspect(self.inputs, 'paired')
        self.assertEqual(paired['read_units'], 100)
        self.assertEqual(paired['read_count'], 200)
        self.assertEqual(rna.inspect(self.inputs[:1], 'single')['read_count'], 100)

    def test_mate_ids_counts_and_orientation_are_validated(self):
        bad = self.base / 'bad.fq'
        for contents, reason in (
            ('@wrong/2\nACGT\n+\nIIII\n', 'IDs disagree'),
            ('@read0/2\nACGT\n+\nIIII\n', 'different read counts'),
            ('@read0/1\nACGT\n+\nIIII\n', 'suffix disagrees')):
            bad.write_text(contents)
            with self.assertRaisesRegex(ValueError, reason):
                rna.inspect([self.inputs[0], bad], 'paired')

    def test_invalid_and_empty_fastq(self):
        bad = self.base / 'bad.fq'
        for contents in ('', '@read\nACGT\n+\nII\n', 'read\nACGT\n+\nIIII\n',
                         '@read\nACGT\n-\nIIII\n', '@read\nACG?\n+\nIIII\n'):
            bad.write_text(contents)
            with self.assertRaises(ValueError):
                rna.inspect([bad], 'single')

    def test_disabled_trimming_is_passthrough(self):
        raw = self.raw()
        clean = self.base / 'clean'
        result = rna.trim('R', {'trim_adapters': False, 'trim_quality': False}, raw, clean, 1, self.inputs)
        self.assertEqual(result['retained_fraction'], 1)
        self.assertTrue(result['passthrough'])
        self.assertEqual(list(clean.iterdir()), [clean / 'rna.json'])
        self.assertEqual(list(raw.iterdir()), [raw / 'rna.json'])
        self.assertEqual(result['sha256'], rna.read_json(raw / 'rna.json')['sha256'])

    def test_fastp_switches_are_explicit(self):
        raw = self.raw()
        settings = {'trim_adapters': True, 'trim_quality': True, 'quality': 20, 'window': 4, 'min_length': 30, 'max_n': 5}
        command, files = rna.fastp_command(raw, settings, self.base, 4, self.inputs)
        for flag in ('--detect_adapter_for_pe', '--cut_tail', '--disable_trim_poly_g', '--disable_quality_filtering'):
            self.assertIn(flag, command)
        self.assertNotIn('--dedup', command)
        self.assertEqual(len(files), 2)
        command, _ = rna.fastp_command(raw, {**settings, 'trim_adapters': False, 'trim_quality': False}, self.base, 1, self.inputs)
        self.assertIn('--disable_adapter_trimming', command)
        self.assertNotIn('--cut_tail', command)

    def test_subsetting_is_reproducible_preserves_pairs_and_full_clean_reads(self):
        raw = self.raw()
        settings = {'subset_fraction': 0.5, 'subset_seed': 42}
        first, second = self.base / 'first', self.base / 'second'
        a = rna.subset('R', settings, raw, first, self.inputs)
        b = rna.subset('R', settings, raw, second, self.inputs)
        self.assertEqual(a['status'], 'ready')
        self.assertEqual(a['read_units'], b['read_units'])
        self.assertGreater(a['read_units'], 0)
        self.assertLess(a['read_units'], 100)
        for mate in (1, 2):
            self.assertEqual(list(rna.records(first / f'R{mate}.fastq.gz')),
                             list(rna.records(second / f'R{mate}.fastq.gz')))
        self.assertEqual(rna.inspect(self.inputs, 'paired')['read_units'], 100)
        full = self.base / 'full'
        self.assertEqual(rna.subset('R', {**settings, 'subset_fraction': 1}, raw, full, self.inputs)['read_units'], 100)
        self.assertEqual(list(full.iterdir()), [full / 'rna.json'])

    def test_failed_source_blocks_only_its_consumers_and_auto_strand_remains_pending(self):
        plan = {'readsets': [{'readset_id': rid, 'library_id': rid, 'strandedness': 'auto'} for rid in ('bad', 'good')],
                'assemblies': [{'assembly_id': aid, 'rna_readsets': reads}
                               for aid, reads in (('A', ['bad']), ('B', ['bad']), ('C', ['good']), ('D', []))]}
        result = report.summarize(plan, [{'readset_id': 'bad', 'status': 'failed'},
                                        {'readset_id': 'good', 'status': 'ready', 'files': ['R1.fastq.gz']}], self.base)
        self.assertEqual([row['status'] for row in result['assemblies']], ['blocked', 'blocked', 'ready', 'not_requested'])
        self.assertEqual(result['readsets'][1]['strand_status'], 'pending_reference_alignment')
        self.assertEqual(result['status'], 'incomplete')

    def test_accession_acquisition_uses_same_integrity_route(self):
        def fake_run(command, **kwargs):
            if '--version' in command:
                return subprocess.CompletedProcess(command, 0, stdout='fasterq-dump 3.1.1')
            if command[0] == 'fasterq-dump':
                folder = Path(command[command.index('-O') + 1])
                for mate in (1, 2):
                    (folder / f'SRR123_{mate}.fastq').write_text(f'@run/{mate}\nACGT\n+\nIIII\n')
            return subprocess.CompletedProcess(command, 0)
        with patch.object(rna.subprocess, 'run', side_effect=fake_run):
            result = rna.acquire({**self.request, 'run_accession': 'SRR123'}, [], self.base / 'accession', 2)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['read_units'], 1)
        self.assertTrue((self.base / 'accession/prefetch_command.json').is_file())


if __name__ == '__main__':
    unittest.main()
