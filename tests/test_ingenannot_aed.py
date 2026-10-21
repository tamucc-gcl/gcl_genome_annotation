import csv
import json
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from ingenannot_inputs import check
from final_annotation import aed
from annotation_fixtures import native_aed_fixture
from prepare_protein_evidence import checksum



class NativeAedContracts(unittest.TestCase):
    def test_preparation_requires_matching_final_annotation_and_evidence_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gff,validation,identity,native=native_aed_fixture(root)
            with self.assertRaisesRegex(ValueError,'availability'):
                check(gff,validation,None,root/'protein.gff',root/'missing.json')
            state=json.loads(validation.read_text());state['stage']='braker_prediction'
            validation.write_text(json.dumps(state))
            with self.assertRaisesRegex(ValueError,'validated final annotation'):
                check(gff,validation,root/'rna.gtf',root/'protein.gff',root/'bad_stage.json')
            state['stage']='final_annotation';validation.write_text(json.dumps(state))
            gff.write_text(gff.read_text()+'# changed\n')
            with self.assertRaisesRegex(ValueError,'models differ'):
                check(gff,validation,root/'rna.gtf',root/'protein.gff',root/'changed.json')

    def test_preparation_runs_without_reporting_or_scoring_helpers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
            scripts=root/'helpers'; scripts.mkdir()
            source=Path(__file__).resolve().parents[1]/'py_scripts'
            for name in ('ingenannot_inputs.py','prepare_protein_evidence.py'):
                shutil.copy2(source/name,scripts/name)
            result=subprocess.run([sys.executable,str(scripts/'ingenannot_inputs.py'),
                '--gff',str(gff),'--validation',str(v),
                '--rna',str(root/'rna.gtf'),'--protein',str(root/'protein.gff'),
                '--out',str(root/'independent_inputs.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(json.loads((root/'independent_inputs.json').read_text()),
                             json.loads((root/'inputs.json').read_text()))

    def test_native_scores_are_reported_not_recomputed_and_all_modes_preserved(self):
        for r,p in [(True,True),(True,False),(False,True),(False,False)]:
            with self.subTest(rna=r,protein=p),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); gff,v,i,a=native_aed_fixture(root,r,p)
                data=aed(gff,a,i,root/'inputs.json',True,.25,root/'out')
                self.assertTrue(data['models_unchanged'])
                with (root/'out/aed_scores.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(row['rna_aed'],'0.3' if r else '.')
                self.assertEqual(row['protein_aed'],'0.1' if p else '.')
                self.assertFalse(any('combined' in k for k in row))

    def test_changed_models_sources_and_ids_block_assessment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
            a.write_text(a.read_text().replace('exon\t101\t200','exon\t101\t199'))
            with self.assertRaisesRegex(ValueError,'changed model structure'): aed(gff,a,i,root/'inputs.json',True,.25,root/'out')

    def test_native_score_outside_range_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
            a.write_text(a.read_text().replace('aed_ev_tr=0.3000','aed_ev_tr=1.3000'))
            with self.assertRaisesRegex(ValueError,'Invalid native AED'): aed(gff,a,i,root/'inputs.json',True,.25,root/'out')

    def test_score_one_is_not_reported_as_zero_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
            a.write_text(a.read_text().replace('aed_ev_tr=0.3000','aed_ev_tr=1.0000').replace('ev_tr=r1','ev_tr=None'))
            data=aed(gff,a,i,root/'inputs.json',True,.25,root/'out')
            stats=data['metrics']['rna']
            self.assertNotIn('no_overlap',stats)
            self.assertEqual(stats['score_equals_one'],1)
            with (root/'out/aed_scores.tsv').open() as h: row=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(row['rna_evidence_id'],'None')

    def test_exon_envelope_normalization_preserves_blocks_on_both_strands(self):
        for strand in ('+','-'):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
                # Mirror BRAKER's separate terminal stop codon: mRNA extends
                # beyond exons/CDS, native export recomputes the exon envelope.
                def rewrite(path,native):
                    lines=[]
                    for line in path.read_text().splitlines():
                        fields=line.split('\t')
                        if len(fields)==9:
                            fields[6]=strand
                            if fields[2] in ('gene','mRNA','transcript'):
                                fields[3]='101' if native and fields[2]!='gene' else '98'
                            line='\t'.join(fields)
                        lines.append(line)
                    path.write_text('\n'.join(lines)+'\n')
                rewrite(gff,False); rewrite(a,True)
                accepted=json.loads((root/'inputs.json').read_text()); accepted['model_sha256']=checksum(gff)
                (root/'inputs.json').write_text(json.dumps(accepted))
                data=aed(gff,a,i,root/'inputs.json',True,.25,root/'out')
                self.assertEqual(data['transcript_bounds_normalized_to_exons'],['t1'])
                a.write_text(a.read_text().replace('mRNA\t101','mRNA\t99'))
                with self.assertRaisesRegex(ValueError,'unexpected transcript bounds'):
                    aed(gff,a,i,root/'inputs.json',True,.25,root/'out')

    def test_cds_phase_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,i,a=native_aed_fixture(root)
            lines=[]
            for line in a.read_text().splitlines():
                fields=line.split('\t')
                if len(fields)==9 and fields[2]=='CDS':
                    fields[7]='1' if fields[7]!='1' else '2'
                    line='\t'.join(fields)
                lines.append(line)
            a.write_text('\n'.join(lines)+'\n')
            with self.assertRaisesRegex(ValueError,'coding phases differ'):
                aed(gff,a,i,root/'inputs.json',True,.25,root/'out')


if __name__=='__main__': unittest.main()
