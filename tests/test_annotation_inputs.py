"""Future cluster unit tests: python -m unittest discover -s tests -p 'test_annotation_inputs.py'."""
import importlib.util
import json
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
            [tid, parent, rank, '', '0', '0', '1', '0', '2', '0', '0', '0', '']) + '\t|\n'
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
