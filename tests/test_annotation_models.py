import sys
from pathlib import Path
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from annotation_models import models, features, rna_hits, merged, overlap_length
from annotation_fixtures import MODEL

class AnnotationModelContracts(unittest.TestCase):


    def test_interval_unions_and_overlap(self):
        self.assertEqual(merged([(10,20),(15,30),(10,20),(40,50)]),[(10,30),(40,50)])
        self.assertEqual(overlap_length([(10,30)],[(20,40)]),10)
        with self.assertRaises(ValueError):merged([(4,4)])

    def test_reference_strand_parent_and_bounds_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.gff3'
            for text in (MODEL.replace('chr1','missing'), MODEL.replace('exon\t101\t200','exon\t99\t200'),
                         MODEL.replace('Parent=t1','Parent=missing'), MODEL.replace('\t+\t','\t.\t')):
                path.write_text(text)
                with self.assertRaises(ValueError):models(path,{'chr1':1000})

    def test_rna_gtf_bounds_and_ambiguous_strand_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'rna.gtf'
            text='chr1\tStringTie\ttranscript\t101\t400\t.\t.\t.\ttranscript_id "r1";\n'+\
                 'chr1\tStringTie\texon\t101\t200\t.\t.\t.\ttranscript_id "r1";\n'
            path.write_text(text)
            self.assertEqual(list(rna_hits(path,{'chr1':1000}))[0]['strand'],'.')
            path.write_text(text.replace('exon\t101\t200','exon\t99\t200'))
            with self.assertRaisesRegex(ValueError,'outside transcript'):list(rna_hits(path,{'chr1':1000}))

if __name__=='__main__':unittest.main()
