import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from review_utr_candidates import classify, review


class UtrReview(unittest.TestCase):
    def test_classification_distinguishes_records_from_base_changes(self):
        self.assertEqual(classify([(10, 20), (10, 20)], [(10, 20)]), 'duplicate_exon_records_only')
        self.assertEqual(classify([(10, 15), (15, 20)], [(10, 20)]), 'same_exonic_bases_different_segmentation')
        self.assertEqual(classify([(10, 20)], [(5, 20)]), 'added_exonic_bases_only')
        self.assertEqual(classify([(10, 20)], [(15, 20)]), 'removed_exonic_bases_only')
        self.assertEqual(classify([(10, 20)], [(5, 15)]), 'added_and_removed_exonic_bases')

    def test_review_reports_restoration_and_selected_rna_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); base = root/'structural_annotation'
            original = base/'braker/A/A.gff3'; final = base/'utr/A/A.utr.gff3'
            native = final.parent/'native/native.gff3'
            identity = root/'filtered_assemblies/A/A.contigs.json'
            rna = original.parent/'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
            for path in (original, final, native, identity, rna): path.parent.mkdir(parents=True, exist_ok=True)
            def gff(extended=False, duplicate=False):
                rows = []
                for n, lo, hi in [('1', 11, 20), ('2', 31, 40)]:
                    start = 6 if extended and n=='1' else lo
                    tag = ';utr_refine_evidence=r1' if extended and n=='1' else ''
                    rows.extend([f'chr1\tx\tgene\t{start}\t{hi}\t.\t+\t.\tID=g{n}',
                        f'chr1\tx\tmRNA\t{start}\t{hi}\t.\t+\t.\tID=t{n};Parent=g{n}{tag}',
                        f'chr1\tx\texon\t{start}\t{hi}\t.\t+\t.\tID=e{n};Parent=t{n}',
                        f'chr1\tx\tCDS\t{lo}\t{hi}\t.\t+\t0\tID=c{n};Parent=t{n}'])
                    if duplicate and n=='2':
                        rows.append(f'chr1\tx\texon\t{lo}\t{hi}\t.\t+\t.\tID=e2dup;Parent=t2')
                return '\n'.join(rows)+'\n'
            original.write_text(gff(duplicate=True)); native.write_text(gff(extended=True))
            final.write_text(gff(extended=True, duplicate=True))
            identity.write_text(json.dumps(dict(reference_identity='ref', contigs=[dict(id='chr1',length=100)])))
            state = dict(status='ready',reference_identity='ref',
                original_model_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
                gff3_sha256=hashlib.sha256(final.read_bytes()).hexdigest(),
                rejected_utr_candidates={'t2':'missing provenance'})
            (final.parent/'utr_status.json').write_text(json.dumps(state))
            with (final.parent/'utr_changes.tsv').open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=['transcript_id','changed','five_prime_utr_bp','three_prime_utr_bp'],delimiter='\t')
                writer.writeheader();writer.writerows([
                    dict(transcript_id='t1',changed=True,five_prime_utr_bp=5,three_prime_utr_bp=0),
                    dict(transcript_id='t2',changed=False,five_prime_utr_bp=0,three_prime_utr_bp=0)])
            rna.write_text('chr1\tx\ttranscript\t6\t20\t.\t+\t.\ttranscript_id "r1";\n'
                           'chr1\tx\texon\t6\t20\t.\t+\t.\ttranscript_id "r1";\n')
            result=review(root,'A',root/'review')
            self.assertEqual(result['rejected_categories'],{'duplicate_exon_records_only':1})
            self.assertTrue(result['all_rejected_originals_restored'])
            case=next(c for c in result['selected_cases'] if c['transcript_id']=='t1')
            self.assertEqual(case['extension_bases_outside_selected_rna'],0)
            self.assertEqual(case['accepted_extension'],[(5,10)])
