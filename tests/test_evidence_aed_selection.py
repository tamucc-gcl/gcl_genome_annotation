import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from evidence_aed import normalize
from evidence_aed_selection import score, select_chains
from test_evidence_aed import IDENTITY, MODEL


def hit(identifier, blocks, target=None):
    return dict(identifier=identifier, target=target or identifier, blocks=blocks)


class IntactChainSelection(unittest.TestCase):
    def test_long_decoy_chain_no_longer_dilutes_a_matching_chain(self):
        picked, count = select_chains([(100,200)], [hit('good',[(100,200)]), hit('decoy',[(100,110),(300,800)])])
        self.assertEqual(picked['identifier'], 'good')
        self.assertEqual(picked['metrics']['aed'], 0)
        self.assertEqual(count, 2)

    def test_no_model_clipping_or_fusion_hiding(self):
        picked, _ = select_chains([(100,200)], [hit('fusion',[(50,250)])])
        self.assertEqual(picked['metrics']['aed'], .25)
        self.assertEqual(picked['blocks'], ((50,250),))

    def test_incompatible_partial_isoforms_are_not_invented_into_one_chain(self):
        picked, count = select_chains([(100,200),(300,400)],
            [hit('left',[(100,200)]),hit('right',[(300,400)])])
        self.assertEqual(picked['metrics']['aed'], .25)
        self.assertEqual(count, 2)

    def test_duplicates_and_iteration_order_do_not_affect_choice(self):
        hits = [hit('a',[(100,200)],'z_source'),hit('b',[(100,200)],'a_source')]
        a, count = select_chains([(100,200)], hits)
        b, reverse_count = select_chains([(100,200)], list(reversed(hits)))
        self.assertEqual(a,b)
        self.assertEqual(a['target'],'a_source')
        self.assertEqual(count,reverse_count)
        self.assertEqual(count,1)

    def test_intronic_span_match_has_no_selected_support(self):
        picked, count = select_chains([(100,200),(300,400)], [hit('intronic',[(210,290)])])
        self.assertIsNone(picked)
        self.assertEqual(count,0)

    def test_selected_conflicting_overhangs_still_penalize_combined_agreement(self):
        from evidence_aed import nucleotide_aed
        rna, _ = select_chains([(100,200)],[hit('rna',[(80,200)])])
        protein, _ = select_chains([(100,200)],[hit('protein',[(100,240)])])
        result = nucleotide_aed([(100,200)],list(rna['blocks'])+list(protein['blocks']))
        self.assertGreater(result['aed'],rna['metrics']['aed'])
        self.assertGreater(result['aed'],protein['metrics']['aed'])

    def test_all_modes_complementary_sources_and_original_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root/'models.gff3'; gff.write_text(MODEL)
            rna = root/'rna.gtf'; rna.write_text(
                'chr1\tStringTie\ttranscript\t101\t200\t.\t+\t.\ttranscript_id "r1";\n'
                'chr1\tStringTie\texon\t101\t200\t.\t+\t.\ttranscript_id "r1";\n')
            protein = root/'protein.gff'; protein.write_text(
                'chr1\tminiprot\tmRNA\t301\t400\t.\t+\t.\tID=p1;Target=source_protein 1 34\n'
                'chr1\tminiprot\tCDS\t301\t400\t.\t+\t0\tParent=p1\n')
            for name,r,p,expected in [('both',rna,protein,0),('rna',rna,None,.25),('protein',None,protein,.25),('none',None,None,None)]:
                with self.subTest(mode=name):
                    normalize(IDENTITY,r,p,root/name)
                    result = score(gff,IDENTITY,root/name,root/(name+'_out'))
                    self.assertEqual(result['gene_aed']['g1'],expected)
                    with (root/(name+'_out')/'aed_scores.tsv').open() as handle:
                        row = next(csv.DictReader(handle,delimiter='\t'))
                    self.assertEqual(row['aed'],row['combined_aed'])
                    self.assertEqual(row['combined_aed'],row['all_union_combined_aed'])
                    self.assertEqual(row['protein_best_source_id'],'source_protein' if p else '.')

    def test_comparison_preserves_original_union_and_no_schema_change_to_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root/'models.gff3'; gff.write_text(MODEL)
            rna = root/'rna.gtf'; rna.write_text(
                'chr1\tStringTie\ttranscript\t101\t400\t.\t+\t.\ttranscript_id "good";\n'
                'chr1\tStringTie\texon\t101\t200\t.\t+\t.\ttranscript_id "good";\n'
                'chr1\tStringTie\texon\t301\t400\t.\t+\t.\ttranscript_id "good";\n'
                'chr1\tStringTie\ttranscript\t1\t1000\t.\t+\t.\ttranscript_id "decoy";\n'
                'chr1\tStringTie\texon\t1\t1000\t.\t+\t.\ttranscript_id "decoy";\n')
            normalize(IDENTITY,rna,None,root/'db')
            before = (root/'db/evidence.sqlite').read_bytes()
            score(gff,IDENTITY,root/'db',root/'best')
            score(gff,IDENTITY,root/'db',root/'union',selection='all_union')
            with (root/'best/aed_scores.tsv').open() as handle: a = next(csv.DictReader(handle,delimiter='\t'))
            with (root/'union/aed_scores.tsv').open() as handle: b = next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(a['combined_aed'],'0.00000000')
            self.assertEqual(a['all_union_combined_aed'],'0.40000000')
            self.assertEqual(a['all_union_combined_aed'],b['combined_aed'])
            self.assertEqual((root/'db/evidence.sqlite').read_bytes(),before)
            self.assertEqual(gff.read_text(),MODEL)

    def test_validation_and_checksum_failures_still_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root/'models.gff3'; gff.write_text(MODEL)
            normalize(IDENTITY,None,None,root/'db')
            with self.assertRaisesRegex(ValueError,'BRAKER validation'):
                score(gff,IDENTITY,root/'db',root/'out',dict(status='ready',gff3_sha256='wrong',reference_identity='test-ref'))
            wrong = dict(IDENTITY,reference_identity='other')
            with self.assertRaisesRegex(ValueError,'reference/checksum'):
                score(gff,wrong,root/'db',root/'out')


if __name__ == '__main__':
    unittest.main()
