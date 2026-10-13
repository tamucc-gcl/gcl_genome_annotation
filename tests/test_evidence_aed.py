import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from evidence_aed import nucleotide_aed, normalize, score, protein_hits, models
from prepare_protein_evidence import checksum


IDENTITY = {'reference_identity': 'test-ref', 'contigs': [{'id': 'chr1', 'length': 1000}]}
MODEL = '''##gff-version 3
chr1\tBRAKER\tgene\t101\t400\t.\t+\t.\tID=g1
chr1\tBRAKER\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1
chr1\tBRAKER\texon\t101\t200\t.\t+\t.\tParent=t1
chr1\tBRAKER\texon\t301\t400\t.\t+\t.\tParent=t1
chr1\tBRAKER\tCDS\t151\t200\t.\t+\t0\tParent=t1
chr1\tBRAKER\tCDS\t301\t350\t.\t+\t1\tParent=t1
'''


class AEDContracts(unittest.TestCase):
    def test_complementary_sources_are_unioned_not_averaged_or_best_hit(self):
        model = [(100, 200), (300, 400)]
        rna, protein = [(100, 200)], [(300, 400)]
        self.assertEqual(nucleotide_aed(model, rna)['aed'], .25)
        self.assertEqual(nucleotide_aed(model, protein)['aed'], .25)
        self.assertEqual(nucleotide_aed(model, rna + protein)['aed'], 0)

    def test_duplicates_do_not_add_weight_and_utr_extensions_are_not_clipped(self):
        self.assertEqual(nucleotide_aed([(100, 200)], [(100, 200)] * 10)['aed'], 0)
        result = nucleotide_aed([(100, 200)], [(50, 250)])
        self.assertEqual(result['aed'], .25)
        self.assertEqual(result['sensitivity'], .5)
        self.assertEqual(result['specificity'], 1)

    def test_intronic_conflict_and_no_support(self):
        self.assertEqual(nucleotide_aed([(100, 200), (300, 400)], [(200, 300)])['aed'], 1)
        self.assertEqual(nucleotide_aed([(100, 200)], [])['aed'], 1)
        with self.assertRaises(ValueError): nucleotide_aed([], [])

    def test_matches_nucleotide_vector_reference_for_random_overlap_patterns(self):
        # Independent bit-set oracle for pinned MAKER get_SN_SP arithmetic.
        rng = random.Random(42)
        for _ in range(200):
            model = [(s, s + rng.randrange(1, 15)) for s in (rng.randrange(100), rng.randrange(100))]
            evidence = [(s, s + rng.randrange(1, 15)) for s in [rng.randrange(100) for _ in range(6)]]
            a = set().union(*(set(range(s, e)) for s, e in model))
            b = set().union(*(set(range(s, e)) for s, e in evidence))
            expected = 1 - (len(a & b) / len(a) + len(a & b) / len(b)) / 2
            self.assertAlmostEqual(nucleotide_aed(model, evidence)['aed'], expected)

    def test_evidence_modes_and_unchanged_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            models = folder / 'models.gff3'; models.write_text(MODEL)
            rna = folder / 'rna.gtf'; rna.write_text(
                'chr1\tStringTie\ttranscript\t101\t200\t.\t+\t.\ttranscript_id "r1";\n'
                'chr1\tStringTie\texon\t101\t200\t.\t+\t.\ttranscript_id "r1";\n')
            protein = folder / 'protein.gff'; protein.write_text(
                'chr1\tminiprot\tmRNA\t301\t400\t.\t+\t.\tID=p1;Target=reference1 1 34\n'
                'chr1\tminiprot\tCDS\t301\t400\t.\t+\t0\tParent=p1\n')
            original = models.read_bytes()
            for name, r, p, expected in [('both', rna, protein, 0), ('rna', rna, None, .25),
                                          ('protein', None, protein, .25), ('none', None, None, None)]:
                with self.subTest(mode=name):
                    normalized = folder / name
                    normalize(IDENTITY, r, p, normalized)
                    result = score(models, IDENTITY, normalized, folder / (name + '_score'))
                    self.assertEqual(result['gene_aed']['g1'], expected)
                    self.assertEqual(models.read_bytes(), original)
                    self.assertEqual(result['transcripts'], 1)

    def test_opposite_strand_evidence_is_not_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            models = folder / 'models.gff3'; models.write_text(MODEL)
            rna = folder / 'rna.gtf'; rna.write_text(
                'chr1\tStringTie\ttranscript\t101\t400\t.\t-\t.\ttranscript_id "r1";\n'
                'chr1\tStringTie\texon\t101\t400\t.\t-\t.\ttranscript_id "r1";\n')
            normalize(IDENTITY, rna, None, folder / 'db')
            result = score(models, IDENTITY, folder / 'db', folder / 'result')
            self.assertEqual(result['gene_aed']['g1'], 1)

    def test_unknown_reference_and_malformed_protein_blocks_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'protein.gff'
            for text in ('chr2\tminiprot\tmRNA\t1\t30\t.\t+\t.\tID=p;Target=q 1 10\n',
                         'chr1\tminiprot\tCDS\t1\t30\t.\t+\t0\tParent=missing\n'):
                path.write_text(text)
                with self.assertRaises(ValueError): list(protein_hits(path, {'chr1': 1000}))

    def test_ambiguous_rna_is_counted_and_empty_supplied_source_is_not_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff = folder / 'models.gff3'; gff.write_text(MODEL)
            rna = folder / 'rna.gtf'; rna.write_text(
                'chr1\tStringTie\ttranscript\t101\t400\t.\t.\t.\ttranscript_id "r1";\n'
                'chr1\tStringTie\texon\t101\t400\t.\t.\t.\ttranscript_id "r1";\n')
            normalized = normalize(IDENTITY, rna, None, folder / 'db')
            self.assertEqual(normalized['evidence']['rna']['ambiguous_strand_excluded'], 1)
            self.assertEqual(score(gff, IDENTITY, folder / 'db', folder / 'out')['gene_aed']['g1'], 1)

    def test_changed_models_and_normalized_database_are_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff = folder / 'models.gff3'; gff.write_text(MODEL)
            normalize(IDENTITY, None, None, folder / 'db')
            accepted = dict(status='ready', gff3_sha256=checksum(gff), reference_identity='test-ref')
            score(gff, IDENTITY, folder / 'db', folder / 'out', accepted)
            gff.write_text(MODEL + '# changed\n')
            with self.assertRaisesRegex(ValueError, 'BRAKER validation'):
                score(gff, IDENTITY, folder / 'db', folder / 'out', accepted)
            with (folder / 'db/evidence.sqlite').open('ab') as handle: handle.write(b'changed')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                score(gff, IDENTITY, folder / 'db', folder / 'out')

    def test_shared_exons_have_explicit_multiple_transcript_parents(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / 'models.gff3'
            gff.write_text(MODEL.replace('Parent=t1\n', 'Parent=t1,t2\n') +
                           'chr1\tBRAKER\tmRNA\t101\t400\t.\t+\t.\tID=t2;Parent=g1\n')
            parsed = models(gff, {'chr1': 1000})
            self.assertEqual(parsed['t1']['exons'], parsed['t2']['exons'])


if __name__ == '__main__':
    unittest.main()
