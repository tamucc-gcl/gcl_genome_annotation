import csv
import json
from pathlib import Path
import tempfile
import unittest
import sys
import shutil
import subprocess

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))

from final_annotation import aed
from final_annotation_inputs import finalize
from ingenannot_inputs import check
from prepare_protein_evidence import checksum
import test_ingenannot_benchmark as benchmark_tests


class FinalAnnotationContracts(unittest.TestCase):
    def test_final_handoff_retains_coding_files_and_blocks_mutated_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gff,v,current,identity,baseline,native=benchmark_tests.NativeBenchmarkContracts().fixtures(root)
            (root/'A.gff3').write_bytes(gff.read_bytes())
            proteins=root/'A_proteins.fa';proteins.write_text('>t1\nMK\n')
            cds=root/'A_cds.fa';cds.write_text('>t1\nATGAAA\n')
            source=dict(status='ready',assembly_id='A',genes=1,coding_transcripts=1,
                gff3_sha256=checksum(gff),proteins_sha256=checksum(proteins),cds_sha256=checksum(cds),reference_identity='test-ref')
            v.write_text(json.dumps(source))
            utr=root/'utr.json';utr.write_text(json.dumps(dict(status='ready',original_model_sha256=checksum(gff),
                gff3_sha256=checksum(gff),reference_identity='test-ref',skipped=False,changed_transcripts=0,
                cds_coordinates_and_phases_unchanged=True,cds_sequences_unchanged=True,protein_sequences_unchanged=True)))
            result=finalize('A',root,v,gff,utr,root/'final',mode='rna_protein')
            self.assertEqual(result['stage'],'final_annotation')
            self.assertEqual((root/'final/A_proteins.fa').read_bytes(),proteins.read_bytes())
            helpers=root/'independent';helpers.mkdir()
            for name in ('final_annotation_inputs.py','prepare_protein_evidence.py'):
                shutil.copy2(Path(__file__).resolve().parents[1]/'py_scripts'/name,helpers/name)
            command=[sys.executable,str(helpers/'final_annotation_inputs.py'),'--assembly-id','A','--mode','rna_protein',
                '--prediction',str(root),'--validation',str(v),'--gff',str(gff),'--utr-status',str(utr),'--outdir',str(root/'independent_result')]
            executed=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(executed.returncode,0,executed.stderr)
            self.assertEqual(check(root/'final/A.gff3',root/'final/annotation_status.json',None,
                root/'rna.gtf',root/'protein.gff',root/'production_inputs.json')['status'],'ready')
            with self.assertRaisesRegex(ValueError,'availability'):
                check(root/'final/A.gff3',root/'final/annotation_status.json',None,None,root/'protein.gff',root/'bad.json')
            proteins.write_text('>t1\nCHANGED\n')
            with self.assertRaisesRegex(ValueError,'coding sequences changed'):
                finalize('A',root,v,gff,utr,root/'bad',mode='rna_protein')

    def test_native_final_scores_handle_all_evidence_modes_without_aliases(self):
        for rna,protein in ((True,True),(True,False),(False,True),(False,False)):
            with self.subTest(rna=rna,protein=protein),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);gff,v,c,identity,baseline,native=benchmark_tests.NativeBenchmarkContracts().fixtures(root,rna,protein)
                state=json.loads(v.read_text());state.update(stage='final_annotation',
                    mode='rna_protein' if rna and protein else 'rna' if rna else 'protein' if protein else 'ab_initio')
                v.write_text(json.dumps(state))
                check(gff,v,None,root/'rna.gtf' if rna else None,root/'protein.gff' if protein else None,root/'inputs.json')
                result=aed(gff,native,identity,root/'inputs.json',True,.25,root/'final')
                with (root/'final/aed_scores.tsv').open() as handle:row=next(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(row['rna_aed'],'0.3' if rna else '.')
                self.assertEqual(row['protein_aed'],'0.1' if protein else '.')
                self.assertNotIn('combined_aed',row)
                self.assertEqual(result['model_sha256'],checksum(gff))
                # Native scoring is assessment only; structure edits remain fatal.
                native.write_text(native.read_text().replace('\tCDS\t151','\tCDS\t152'))
                with self.assertRaisesRegex(ValueError,'changed model structure'):
                    aed(gff,native,identity,root/'inputs.json',True,.25,root/'bad')
