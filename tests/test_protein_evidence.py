"""Cluster unit tests for prediction protein acquisition and source accounting."""
import gzip
import importlib.util
import io
import json
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('proteins', ROOT / 'py_scripts/prepare_protein_evidence.py')
proteins = importlib.util.module_from_spec(spec)
spec.loader.exec_module(proteins)
report_spec = importlib.util.spec_from_file_location('protein_report', ROOT / 'py_scripts/report_protein_evidence.py')
protein_report = importlib.util.module_from_spec(report_spec)
report_spec.loader.exec_module(protein_report)


class ProteinTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.local = self.base / 'spike.fa'
        self.local.write_text('>shared description\nmpeptidex*\n>other\nMKL\n')
        self.catalogue = json.loads((ROOT / 'assets/orthodb_sources.json').read_text())

    def row(self, aid='A', **overrides):
        return {'assembly_id': aid, 'settings': {'use_protein_evidence': True,
            'orthodb_clade': 'none', 'ref_proteome_accession': [],
            'protein_fastas': [str(self.local)], **overrides}}

    def plan(self, rows, name='plan', **options):
        return proteins.plan({'status': 'ready', 'assemblies': rows}, self.catalogue,
                             options, self.base / name)

    def source(self, spec, name='sources'):
        directory = self.base / name
        directory.mkdir(exist_ok=True)
        target = directory / spec['source_id']
        return directory, proteins.acquire(spec, [spec['local_path']] if spec.get('local_path') else [], target)

    def test_shared_sources_are_unique_and_assembly_order_is_irrelevant(self):
        first = self.plan([self.row('A'), self.row('B')])
        second = self.plan([self.row('B'), self.row('A')], name='second')
        self.assertEqual(len(first['sources']), 1)
        self.assertEqual(first['sources'], second['sources'])
        self.assertEqual(first['assemblies'][0]['source_ids'], first['assemblies'][1]['source_ids'])

    def test_disabled_proteins_have_no_sentinel(self):
        plan = self.plan([self.row(use_protein_evidence=False)])
        self.assertEqual(plan['sources'], [])
        target = self.base / 'evidence'
        result = proteins.combine(plan['assemblies'][0], self.base / 'absent', target)
        self.assertEqual(result['status'], 'not_requested')
        self.assertFalse((target / 'proteins.fa').exists())

    def test_normalization_and_original_header_crosswalk(self):
        plan = self.plan([self.row()])
        _, source = self.source(plan['sources'][0])
        self.assertEqual(source['status'], 'ready')
        self.assertEqual(source['protein_count'], 2)
        self.assertEqual(source['terminal_stops_removed'], 1)
        target = self.base / 'sources' / source['source_id']
        self.assertIn('MPEPTIDEX', (target / 'proteins.fa').read_text())
        self.assertNotIn('*', (target / 'proteins.fa').read_text())
        self.assertIn('shared description', (target / 'id_map.tsv').read_text())

    def test_bad_fasta_never_leaves_partial_valid_evidence(self):
        for index, raw in enumerate(('MPEPTIDE\n', '>x\n', '>x\nMP*EP\n',
                                      '>x\nMP-EP\n', '>x\nMP\n>x\nEP\n')):
            with self.subTest(raw=raw):
                target = self.base / f'bad_{index}'
                target.mkdir()
                with self.assertRaises(ValueError):
                    proteins.normalize(io.StringIO(raw), 'p123', target)
                self.assertFalse((target / 'proteins.fa').exists())

    def test_failed_shared_source_blocks_consumers_and_preserves_raw(self):
        self.local.write_text('>bad\nMP*EP\n')
        plan = self.plan([self.row('A'), self.row('B')])
        directory, source = self.source(plan['sources'][0])
        self.assertEqual(source['status'], 'failed')
        self.assertTrue((directory / source['source_id'] / 'raw.fa').is_file())
        for row in plan['assemblies']:
            result = proteins.combine(row, directory, self.base / row['assembly_id'])
            self.assertEqual(result['status'], 'blocked')
            self.assertFalse((self.base / row['assembly_id'] / 'proteins.fa').exists())

    def test_combined_sources_keep_duplicate_original_ids_separate(self):
        other = self.base / 'other.fa'
        other.write_text('>shared other source\nMKK\n')
        plan = self.plan([self.row(protein_fastas=[str(self.local), str(other)])])
        directory = None
        for source in plan['sources']:
            directory, status = self.source(source)
            self.assertEqual(status['status'], 'ready')
        result = proteins.combine(plan['assemblies'][0], directory, self.base / 'evidence')
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['protein_count'], 3)
        ids = [header for header, _ in proteins.records(io.StringIO((self.base / 'evidence/proteins.fa').read_text()))]
        self.assertEqual(len(set(ids)), 3)

    def test_changed_local_input_is_rejected(self):
        plan = self.plan([self.row()])
        self.local.write_text('>changed\nMKK\n')
        _, result = self.source(plan['sources'][0])
        self.assertEqual(result['status'], 'failed')
        self.assertIn('changed after planning', result['reason'])

    def test_changed_normalized_output_blocks_combination(self):
        plan = self.plan([self.row()])
        directory, source = self.source(plan['sources'][0])
        (directory / source['source_id'] / 'proteins.fa').write_text('>tampered\nMK\n')
        result = proteins.combine(plan['assemblies'][0], directory, self.base / 'evidence')
        self.assertEqual(result['status'], 'blocked')
        self.assertIn('checksum mismatch', result['reason'])

    def test_cached_orthodb_requires_provider_checksum(self):
        archive = self.base / 'Vertebrata.fa.gz'
        with gzip.open(archive, 'wt') as handle:
            handle.write('>odb1\nMPEPTIDE\n')
        self.catalogue['partitions']['Vertebrata'] = proteins.checksum(archive, 'md5')
        plan = self.plan([self.row(orthodb_clade='Vertebrata', protein_fastas=[])],
                         orthodb_database=str(self.base), launch_dir=str(self.base))
        _, result = self.source(plan['sources'][0])
        self.assertEqual(result['status'], 'ready')
        bad_spec = {**plan['sources'][0], 'expected_md5': '0' * 32}
        _, bad = self.source(bad_spec, name='bad_cache')
        self.assertEqual(bad['status'], 'failed')
        self.assertIn('MD5', bad['reason'])

    def test_ncbi_exact_accession_package_and_missing_proteins(self):
        plan = self.plan([self.row(protein_fastas=[], ref_proteome_accession=['GCF_024713975.1'])])
        source = plan['sources'][0]

        def fake_run(command, **kwargs):
            if command[1] == 'version':
                return subprocess.CompletedProcess(command, 0, stdout='datasets version: 18.35.0\n')
            destination = Path(command[command.index('--filename') + 1])
            with zipfile.ZipFile(destination, 'w') as archive:
                archive.writestr('ncbi_dataset/data/GCF_024713975.1/protein.faa', '>spike\nMPEPTIDE\n')
                archive.writestr('ncbi_dataset/data/assembly_data_report.jsonl', '{"accession":"GCF_024713975.1"}\n')
            return subprocess.CompletedProcess(command, 0)

        with patch.object(proteins.subprocess, 'run', side_effect=fake_run):
            _, result = self.source(source)
        self.assertEqual(result['status'], 'ready')
        self.assertTrue((self.base / 'sources' / source['source_id'] / 'raw.zip').is_file())

        def wrong_run(command, **kwargs):
            if command[1] == 'version':
                return subprocess.CompletedProcess(command, 0, stdout='18.35.0')
            with zipfile.ZipFile(command[command.index('--filename') + 1], 'w') as archive:
                archive.writestr('ncbi_dataset/data/GCF_999999999.1/protein.faa', '>wrong\nMK\n')
            return subprocess.CompletedProcess(command, 0)

        with patch.object(proteins.subprocess, 'run', side_effect=wrong_run):
            _, result = self.source(source, name='wrong_accession')
        self.assertEqual(result['status'], 'failed')
        self.assertIn('exact-accession', result['reason'])

    def test_only_failed_sources_get_new_resume_identity(self):
        first = self.plan([self.row('A'), self.row('B')])
        directory, source = self.source(first['sources'][0])
        ready = self.plan([self.row()], name='ready_plan', published_sources=str(directory))
        self.assertNotIn('retry_nonce', ready['sources'][0])
        source['status'] = 'failed'
        proteins.write_json(directory / source['source_id'] / 'source.json', source)
        retry = self.plan([self.row('A'), self.row('B')], name='retry_plan', published_sources=str(directory))
        self.assertEqual(retry['sources'][0]['source_id'], first['sources'][0]['source_id'])
        self.assertIn('retry_nonce', retry['sources'][0])

    def test_shared_store_survives_new_run_and_detects_tampering(self):
        cache = self.base / 'general_database'
        cache.mkdir()
        archive = cache / 'Vertebrata.fa.gz'
        with gzip.open(archive, 'wt') as handle:
            handle.write('>ortho\nMPEPTIDE\n')
        self.catalogue['partitions']['Vertebrata'] = proteins.checksum(archive, 'md5')
        options = {'orthodb_database': str(cache), 'launch_dir': str(self.base)}
        row = self.row(protein_fastas=[], orthodb_clade='Vertebrata')
        first = self.plan([row], **options)
        _, source = self.source(first['sources'][0])
        bundle = self.base / 'sources' / source['source_id']
        store = Path(first['sources'][0]['store_dir']) / source['source_id']
        store.parent.mkdir(parents=True)
        import shutil
        shutil.copytree(bundle, store)
        proteins.write_json(store / 'STORE_READY.json', {'status': 'ready'})
        archive.unlink()  # Reuse cannot depend on the original archive shortcut.
        second = self.plan([row, self.row('B', protein_fastas=[], orthodb_clade='Vertebrata')],
                           name='different_project', **options)
        self.assertEqual(len(second['sources']), 1)
        self.assertEqual(second['sources'][0]['source_id'], source['source_id'])
        self.assertEqual(second['sources'][0]['stored_path'], str(store.resolve()))
        (store / 'proteins.fa').write_text('>tampered\nMK\n')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            self.plan([row], name='tampered_store', **options)

    def test_completed_run_bundle_is_imported_without_redownload(self):
        options = {'protein_source_store': str(self.base / 'project_store'), 'launch_dir': str(self.base)}
        first = self.plan([self.row()], **options)
        directory, result = self.source(first['sources'][0])
        second = self.plan([self.row()], name='import_plan', published_sources=str(directory), **options)
        self.assertEqual(second['sources'][0]['stored_path'], str((directory / result['source_id']).resolve()))
        self.assertFalse(Path(second['sources'][0]['store_dir']).exists())

    def test_summary_accounts_for_disabled_and_blocked(self):
        plan = self.plan([self.row('disabled', use_protein_evidence=False), self.row('blocked')])
        directories = []
        for row in plan['assemblies']:
            out = self.base / row['assembly_id']
            proteins.combine(row, self.base / 'missing_sources', out)
            directories.append(out)
        result = protein_report.summarize(directories, self.base / 'status.json', self.base / 'README.md')
        self.assertEqual(result['status'], 'incomplete')
        self.assertEqual([x['status'] for x in result['assemblies']], ['blocked', 'not_requested'])


if __name__ == '__main__':
    unittest.main()
