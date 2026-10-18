import json
from pathlib import Path
import tempfile
import unittest

from inspect_native_aed_models import inspect, difference
from test_evidence_aed import MODEL, IDENTITY
from prepare_protein_evidence import checksum


class NativeAedDiagnostic(unittest.TestCase):
    def test_distinguishes_duplicate_segmentation_and_true_base_changes(self):
        self.assertEqual(difference([(10,20)],[(10,20),(10,20)]),'duplicate_records_only')
        self.assertEqual(difference([(10,20)],[(10,15),(15,20)]),'same_bases_different_segmentation')
        self.assertEqual(difference([(10,20)],[(11,20)]),'changed_bases')

    def test_reports_actual_coordinate_delta_without_relaxing_assessment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original=root/'original.gff3';original.write_text(MODEL)
            native=root/'native.gff3';native.write_text(MODEL.replace('exon\t101','exon\t102'))
            identity=root/'identity.json';identity.write_text(json.dumps(IDENTITY))
            validation=root/'validation.json';validation.write_text(json.dumps(dict(status='ready',
                gff3_sha256=checksum(original),reference_identity=IDENTITY['reference_identity'])))
            result=inspect(original,native,identity,validation,root/'diagnostic.json',['t1'])
            self.assertEqual(result['categories'],{'exons:changed_bases':1})
            self.assertEqual(result['examples'][0]['removed_exonic_bases'],[(100,101)])
            self.assertEqual(result['examples'][0]['added_exonic_bases'],[])
