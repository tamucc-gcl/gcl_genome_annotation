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
from annotation_fixtures import native_aed_fixture


class FinalAnnotationContracts(unittest.TestCase):
    def test_duplicate_native_exons_are_audited_without_changing_annotation(self):
        for strand in ('+','-'):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                gff,v,identity,native=native_aed_fixture(root)
                for path in (gff,native):
                    path.write_text(path.read_text().replace('\t+\t','\t'+strand+'\t'))
                accepted=json.loads((root/'inputs.json').read_text())
                accepted['model_sha256']=checksum(gff)
                (root/'inputs.json').write_text(json.dumps(accepted))
                original=gff.read_bytes()
                exon=next(line for line in native.read_text().splitlines() if '\texon\t' in line)
                native.write_text(native.read_text()+exon+'\n'+exon+'\n')
                result=aed(gff,native,identity,root/'inputs.json',True,.25,root/'out')
                audit=result['native_exon_record_normalizations']['t1']
                self.assertEqual(audit['native_records'],audit['original_records']+2)
                self.assertEqual(gff.read_bytes(),original)
                with (root/'out/aed_scores.tsv').open() as handle:
                    self.assertEqual(next(csv.DictReader(handle,delimiter='\t'))['protein_aed'],'0.1')
                self.assertEqual(json.loads((root/'out/aed_status.json').read_text())['native_exon_record_normalizations']['t1'],audit)

    def test_identical_exon_coverage_with_changed_segmentation_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            gff,v,identity,native=native_aed_fixture(root)
            exon=next(line for line in native.read_text().splitlines() if '\texon\t' in line)
            fields=exon.split('\t');middle=(int(fields[3])+int(fields[4]))//2
            left=fields.copy();left[4]=str(middle)
            right=fields.copy();right[3]=str(middle+1);right[8]=right[8].replace('ID=','ID=split_')
            native.write_text(native.read_text().replace(exon,'\t'.join(left)+'\n'+'\t'.join(right)))
            with self.assertRaisesRegex(ValueError,'changed model structure'):
                aed(gff,native,identity,root/'inputs.json',True,.25,root/'bad')

    def test_duplicate_cds_records_remain_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            gff,v,identity,native=native_aed_fixture(root)
            cds=next(line for line in native.read_text().splitlines() if '\tCDS\t' in line)
            native.write_text(native.read_text()+cds+'\n')
            with self.assertRaisesRegex(ValueError,'changed model structure'):
                aed(gff,native,identity,root/'inputs.json',True,.25,root/'bad')

    def test_final_handoff_retains_coding_files_and_blocks_mutated_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gff,v,identity,native=native_aed_fixture(root)
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
            self.assertEqual(check(root/'final/A.gff3',root/'final/annotation_status.json',
                root/'rna.gtf',root/'protein.gff',root/'production_inputs.json')['status'],'ready')
            with self.assertRaisesRegex(ValueError,'availability'):
                check(root/'final/A.gff3',root/'final/annotation_status.json',None,root/'protein.gff',root/'bad.json')
            proteins.write_text('>t1\nCHANGED\n')
            with self.assertRaisesRegex(ValueError,'coding sequences changed'):
                finalize('A',root,v,gff,utr,root/'bad',mode='rna_protein')

    def test_native_final_scores_handle_all_evidence_modes_without_aliases(self):
        for rna,protein in ((True,True),(True,False),(False,True),(False,False)):
            with self.subTest(rna=rna,protein=protein),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);gff,v,identity,native=native_aed_fixture(root,rna,protein)
                state=json.loads(v.read_text());state.update(stage='final_annotation',
                    mode='rna_protein' if rna and protein else 'rna' if rna else 'protein' if protein else 'ab_initio')
                v.write_text(json.dumps(state))
                check(gff,v,root/'rna.gtf' if rna else None,root/'protein.gff' if protein else None,root/'inputs.json')
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
