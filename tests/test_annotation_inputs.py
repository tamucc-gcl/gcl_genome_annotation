"""Future cluster unit tests: python -m unittest discover -s tests -p 'test_annotation_inputs.py'."""
import importlib.util
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('resolver', ROOT / 'py_scripts/resolve_annotation_inputs.py')
resolver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resolver)


class InputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        (self.base / 'genome.fa').write_text('>chr1\nACGT\n')
        (self.base / 'reads.fq').write_text('@r\nAC\n+\nII\n')
        self.rules = json.loads((ROOT / 'assets/taxonomy_rules.json').read_text())
        self.params = {'use_protein_evidence': False, 'genetic_code': 1,
                       'taxa_name': '2759', 'run_busco': False,
                       'run_telomere': False, 'run_omark': False}

    def assembly(self, **extra):
        row = {'assembly_id': 'A', 'assembly_fasta': 'genome.fa', **extra}
        return resolver.resolve_assembly(row, self.params, None, self.rules, self.base, self.base)

    def test_false_zero_override(self):
        self.params['min_contig_length'] = 1000
        result = self.assembly(min_contig_length='0', run_eggnog='false')
        self.assertEqual(result['settings']['min_contig_length'], 0)
        self.assertFalse(result['settings']['run_eggnog'])
        self.assertEqual(result['origins']['min_contig_length'], 'assembly_row')

    def test_all_four_routes(self):
        self.assertEqual(self.assembly()['prediction_mode'], 'ab_initio')
        self.assertEqual(self.assembly(rna_readsets='R')['prediction_mode'], 'rna')
        self.params.update(use_protein_evidence=True, orthodb_clade='Vertebrata')
        self.assertEqual(self.assembly()['prediction_mode'], 'protein')
        self.assertEqual(self.assembly(rna_readsets='R')['prediction_mode'], 'rna_protein')

    def test_alias_conflict(self):
        with self.assertRaises(ValueError):
            self.assembly(assembly_path='different.fa')

    def test_quoted_csv(self):
        path = self.base / 'sheet.csv'
        path.write_text('assembly_fasta,assembly_id\n"path,with,commas.fa",A\n')
        rows = resolver.read_table(path, resolver.ASSEMBLY_FIELDS, ['assembly_id'])
        self.assertEqual(rows[0]['assembly_fasta'], 'path,with,commas.fa')

    def test_bad_headers(self):
        path = self.base / 'sheet.csv'
        path.write_text('assembly_id,assembly_id\nA,B\n')
        with self.assertRaises(ValueError):
            resolver.read_table(path, resolver.ASSEMBLY_FIELDS, ['assembly_id'])

    def table(self, name, header, rows):
        path = self.base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)
        return str(path)

    def manifest(self, assembly_sheet, rna_sheet=None, **extra):
        return resolver.build_manifest({'launch_dir': str(self.base),
            'sample_sheet': assembly_sheet, 'rna_sheet': rna_sheet,
            'params': self.params, **extra}, self.rules)

    def test_malformed_and_unknown_columns(self):
        for content, message in (
            ('assembly_id,assembly_fasta\nA,genome.fa,\n', 'expected 2 fields, found 3'),
            ('assembly_id,assembly_fasta\nA\n', 'expected 2 fields, found 1'),
            ('assembly_id,assemblly_fasta\nA,genome.fa\n', 'unknown columns'),
            ('assembly_fasta\ngenome.fa\n', 'required columns'),
        ):
            with self.subTest(content=content):
                path = self.base / 'bad.csv'
                path.write_text(content)
                with self.assertRaisesRegex(ValueError, message):
                    resolver.read_table(path, resolver.ASSEMBLY_FIELDS, ['assembly_id'])

    def test_missing_reference_header_is_global_error(self):
        sheet = self.table('a.csv', ['assembly_id'], [['A']])
        with self.assertRaisesRegex(ValueError, 'required column assembly_fasta'):
            self.manifest(sheet)

    def test_duplicate_ids_are_global_errors(self):
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta'],
                           [['A', 'genome.fa'], ['A', 'genome.fa']])
        with self.assertRaisesRegex(ValueError, 'Duplicate assembly_id'):
            self.manifest(sheet)
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta'], [['A', 'genome.fa']])
        rna = self.table('r.csv', ['readset_id', 'library_id', 'layout', 'fastq_1'],
                         [['R', 'L', 'single', 'reads.fq'], ['R', 'L', 'single', 'reads.fq']])
        with self.assertRaisesRegex(ValueError, 'Duplicate readset_id'):
            self.manifest(sheet, rna)

    def test_four_routes_in_one_batch_and_multiple_shared_readsets(self):
        self.params.update(use_protein_evidence=True, orthodb_clade='Vertebrata')
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta', 'rna_readsets',
            'use_protein_evidence'], [['combined', 'genome.fa', 'R2;R1', 'true'],
            ['rna', 'genome.fa', 'R1', 'false'], ['protein', 'genome.fa', '', 'true'],
            ['genome', 'genome.fa', '', 'false']])
        rna = self.table('r.csv', ['readset_id', 'library_id', 'layout', 'fastq_1', 'run_accession'],
            [['R1', 'L1', 'single', 'reads.fq', ''], ['R2', 'L2', 'paired', '', 'ERR123']])
        result = self.manifest(sheet, rna)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual([x['prediction_mode'] for x in result['assemblies']],
                         ['rna_protein', 'rna', 'protein', 'ab_initio'])
        self.assertEqual(result['assemblies'][0]['rna_readsets'], ['R1', 'R2'])
        self.assertEqual(result['assemblies'][2]['rna_readsets'], [])
        self.assertEqual(len(result['rna_readsets']), 2)

    def test_failed_shared_rna_blocks_only_consumers(self):
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta', 'rna_readsets'],
            [['A', 'genome.fa', 'R'], ['B', 'genome.fa', 'R'], ['C', 'genome.fa', '']])
        rna = self.table('r.csv', ['readset_id', 'library_id', 'layout', 'fastq_1'],
                         [['R', 'L', 'single', 'missing.fq']])
        result = self.manifest(sheet, rna)
        self.assertEqual([x['status'] for x in result['assemblies']], ['blocked', 'blocked', 'ready'])
        self.assertEqual(result['rna_readsets'][0]['status'], 'failed')
        self.assertEqual(result['status'], 'incomplete')

    def test_unknown_rna_is_blocked_not_absent(self):
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta', 'rna_readsets'],
                           [['A', 'genome.fa', 'unknown']])
        result = self.manifest(sheet)
        self.assertEqual(result['assemblies'][0]['status'], 'blocked')
        self.assertIn('unknown', result['assemblies'][0]['reason'])

    def test_duplicate_rna_sources_and_inconsistent_libraries(self):
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta'], [['A', 'genome.fa']])
        rna = self.table('r.csv', ['readset_id', 'library_id', 'layout', 'fastq_1'],
            [['R1', 'L', 'single', 'reads.fq'], ['R2', 'L', 'single', 'reads.fq']])
        with self.assertRaisesRegex(ValueError, 'RNA source assigned twice'):
            self.manifest(sheet, rna)
        rna = self.table('r.csv', ['readset_id', 'library_id', 'layout', 'run_accession', 'strandedness'],
            [['R1', 'L', 'single', 'SRR123', 'forward'], ['R2', 'L', 'single', 'SRR124', 'reverse']])
        with self.assertRaisesRegex(ValueError, 'Inconsistent settings within RNA library'):
            self.manifest(sheet, rna)

    def test_paired_single_and_accession_settings(self):
        (self.base / 'mate.fq').write_text('@r\nGT\n+\nII\n')
        row = {'readset_id': 'R', 'library_id': 'L', 'layout': 'paired',
               'fastq_1': 'reads.fq', 'fastq_2': 'mate.fq', 'trim_quality': 'false',
               'subset_fraction': '0.5', 'subset_seed': '0', 'strandedness': 'reverse'}
        result = resolver.resolve_rna(row, {}, self.base)
        self.assertFalse(result['trim_quality'])
        self.assertEqual(result['subset_seed'], 0)
        self.assertEqual(result['subset_fraction'], 0.5)
        for accession in ('SRR123', 'ERR123', 'DRR123'):
            result = resolver.resolve_rna({'readset_id': 'R', 'library_id': 'L',
                'layout': 'single', 'run_accession': accession}, {}, self.base)
            self.assertEqual(result['run_accession'], accession)
        for changes in ({'fastq_2': ''}, {'layout': 'single'}, {'strandedness': 'guess'},
                        {'subset_fraction': '0'}, {'subset_fraction': '1.1'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                resolver.resolve_rna({**row, **changes}, {}, self.base)

    def test_relative_paths_and_legacy_reference_header(self):
        folder = self.base / 'sheets'
        folder.mkdir()
        (folder / 'genome, with spaces.fa').write_text('>chr1\nACGT\n')
        sheet = self.table('sheets/a.csv', ['assembly_path', 'assembly_id'],
                           [['genome, with spaces.fa', 'A']])
        result = self.manifest(sheet)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['assemblies'][0]['assembly_fasta'],
                         str((folder / 'genome, with spaces.fa').resolve()))

    def test_protein_source_inheritance_clear_and_disable(self):
        self.params.update(use_protein_evidence=True, orthodb_clade='none',
                           ref_proteome_accession='GCF_024713975.1', protein_fastas='genome.fa')
        inherited = self.assembly()
        self.assertEqual(inherited['settings']['ref_proteome_accession'], ['GCF_024713975.1'])
        cleared = self.assembly(protein_fastas='none')
        self.assertEqual(cleared['settings']['protein_fastas'], [])
        self.assertEqual(cleared['prediction_mode'], 'protein')
        disabled = self.assembly(use_protein_evidence='false')
        self.assertEqual(disabled['prediction_mode'], 'ab_initio')
        self.assertEqual(disabled['settings']['ref_proteome_accession'], [])
        self.assertEqual(disabled['origins']['ref_proteome_accession'], 'disabled_by_use_protein_evidence')
        self.assertEqual(disabled['origins']['orthodb_clade'], 'disabled_by_use_protein_evidence')
        with self.assertRaisesRegex(ValueError, 'every source is disabled'):
            self.assembly(protein_fastas='none', ref_proteome_accession='none')

    def test_cli_failure_report_and_exit_status(self):
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta'], [['A', 'missing.fa']])
        request = self.base / 'request.json'
        request.write_text(json.dumps({'launch_dir': str(self.base), 'sample_sheet': sheet,
            'rna_sheet': None, 'params': self.params}))
        outdir = self.base / 'report'
        command = [sys.executable, str(ROOT / 'py_scripts/resolve_annotation_inputs.py'),
            '--request', str(request), '--rules', str(ROOT / 'assets/taxonomy_rules.json'),
            '--outdir', str(outdir)]
        failed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(failed.returncode, 2, failed.stderr)
        self.assertEqual(json.loads((outdir / 'resolved_inputs.json').read_text())['status'], 'incomplete')
        self.assertIn('missing.fa', (outdir / 'README.md').read_text())
        report_only = subprocess.run(command + ['--report-only'], capture_output=True, text=True)
        self.assertEqual(report_only.returncode, 0, report_only.stderr)

    def test_taxid_without_snapshot_blocks_assembly(self):
        self.params['taxid'] = '999'
        sheet = self.table('a.csv', ['assembly_id', 'assembly_fasta'], [['A', 'genome.fa']])
        result = self.manifest(sheet)
        self.assertEqual(result['assemblies'][0]['status'], 'blocked')
        self.assertIn('no local --taxdump_dir', result['assemblies'][0]['reason'])

    def test_shared_rna_and_blocked_branch(self):
        (self.base / 'a.csv').write_text('assembly_id,assembly_fasta,rna_readsets\nA,genome.fa,R\nB,genome.fa,R\nC,missing.fa,R\n')
        (self.base / 'r.csv').write_text('readset_id,library_id,layout,fastq_1\nR,L,single,reads.fq\n')
        result = resolver.build_manifest({'launch_dir': str(self.base), 'sample_sheet': 'a.csv',
            'rna_sheet': 'r.csv', 'params': self.params}, self.rules)
        self.assertEqual(len(result['rna_readsets']), 1)
        self.assertEqual([x['status'] for x in result['assemblies']], ['ready', 'ready', 'blocked'])
        self.assertEqual(result['status'], 'incomplete')

    def test_accession_exclusive(self):
        row = {'readset_id': 'R', 'library_id': 'L', 'layout': 'paired',
               'run_accession': 'SRR123', 'fastq_1': 'reads.fq'}
        with self.assertRaises(ValueError):
            resolver.resolve_rna(row, {}, self.base)

    def test_pair_integrity(self):
        with self.assertRaises(ValueError):
            resolver.resolve_rna({'readset_id': 'R', 'library_id': 'L', 'layout': 'paired',
                                  'fastq_1': 'reads.fq', 'fastq_2': 'reads.fq'}, {}, self.base)

    def test_source_conflict(self):
        with self.assertRaises(ValueError):
            self.assembly(protein_fastas='genome.fa')

    def test_no_orthodb_fallback(self):
        self.assertIsNone(resolver.derive('orthodb_clade', {'ancestors': ['2759', '1']}, self.rules))

    def test_global_legacy_taxid(self):
        self.assertEqual(resolver.alias({'ncbi_speciesID': 999}, 'taxid', 'ncbi_speciesID')['taxid'], 999)
        with self.assertRaises(ValueError):
            resolver.alias({'taxid': 999, 'ncbi_speciesID': 1000}, 'taxid', 'ncbi_speciesID')

    def test_taxonomy_precedence_and_nuclear_code(self):
        directory = self.base / 'taxdump'
        directory.mkdir()
        # Synthetic hierarchy includes different mitochondrial and nuclear codes.
        nodes = [('1', '1', 'no rank'), ('2759', '1', 'superkingdom'),
                 ('7742', '2759', 'clade'), ('7898', '7742', 'class'), ('999', '7898', 'species')]
        (directory / 'nodes.dmp').write_text(''.join('\t|\t'.join(
            [tid, parent, rank, '', '0', '0', '6' if tid == '999' else '1',
             '1' if tid == '999' else '0', '2', '0', '0', '0', '']) + '\t|\n'
            for tid, parent, rank in nodes))
        (directory / 'names.dmp').write_text(''.join(f'{tid}\t|\tTaxon {tid}\t|\t\t|\tscientific name\t|\n' for tid, _, _ in nodes))
        (directory / 'merged.dmp').write_text('998\t|\t999\t|\n')
        tax = resolver.Taxonomy(directory)
        self.assertEqual(tax.resolve('998')['taxid'], '999')
        self.params.update(taxid='999', use_protein_evidence=True, genetic_code='auto', taxa_name='auto')
        row = {'assembly_id': 'A', 'assembly_fasta': 'genome.fa'}
        result = resolver.resolve_assembly(row, self.params, tax, self.rules, self.base, self.base)
        self.assertEqual(result['settings']['genetic_code'], 1)
        self.assertEqual(result['settings']['orthodb_clade'], 'Vertebrata')
        self.assertEqual(result['settings']['busco_lineage'], 'actinopterygii_odb10')
        self.params['busco_lineage'] = 'global_override'
        result = resolver.resolve_assembly(row, self.params, tax, self.rules, self.base, self.base)
        self.assertEqual(result['settings']['busco_lineage'], 'global_override')
        row.update(busco_lineage='row_override', taxid='7898')
        result = resolver.resolve_assembly(row, self.params, tax, self.rules, self.base, self.base)
        self.assertEqual(result['settings']['taxid'], '7898')
        self.assertEqual(result['settings']['busco_lineage'], 'row_override')
        row['busco_lineage'] = 'auto'
        result = resolver.resolve_assembly(row, self.params, tax, self.rules, self.base, self.base)
        self.assertEqual(result['settings']['busco_lineage'], 'actinopterygii_odb10')


if __name__ == '__main__':
    unittest.main()
