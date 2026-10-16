import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from ingenannot_benchmark import check, report
from prepare_protein_evidence import checksum
from test_evidence_aed import IDENTITY, MODEL


class NativeBenchmarkContracts(unittest.TestCase):
    def fixtures(self,root,rna=True,protein=True):
        gff=root/'models.gff3'; gff.write_text(MODEL)
        validation=root/'validation.json'; validation.write_text(json.dumps(dict(status='ready',gff3_sha256=checksum(gff),reference_identity='test-ref')))
        r=root/'rna.gtf'; r.write_text('# RNA fixture\n')
        p=root/'protein.gff'; p.write_text('# protein fixture\n')
        current=root/'current.json'; current.write_text(json.dumps(dict(status='ready',model_sha256=checksum(gff),reference_identity='test-ref',
            evidence=dict(rna=dict(available=rna,source_sha256=checksum(r) if rna else None),protein=dict(available=protein,source_sha256=checksum(p) if protein else None)))))
        identity=root/'identity.json'; identity.write_text(json.dumps(IDENTITY))
        baseline=root/'scores.tsv'; baseline.write_text('gene_id\ttranscript_id\trna_aed\tprotein_aed\n'+f'g1\tt1\t{0.4 if rna else "."}\t{0.2 if protein else "."}\n')
        annotated=root/'native.gff3'; annotated.write_text(MODEL.replace('ID=t1;Parent=g1','ID=t1;Parent=g1;aed_ev_tr=0.3000;aed_ev_pr=0.1000;ev_tr=r1;ev_pr=p1;ev_tr_penalty=no'))
        check(gff,validation,current,r if rna else None,p if protein else None,root/'inputs.json')
        return gff,validation,current,identity,baseline,annotated

    def test_native_scores_are_reported_not_recomputed_and_all_modes_preserved(self):
        for r,p in [(True,True),(True,False),(False,True),(False,False)]:
            with self.subTest(rna=r,protein=p),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); gff,v,c,i,b,a=self.fixtures(root,r,p)
                data=report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')
                self.assertTrue(data['models_unchanged'])
                with (root/'out/benchmark_scores.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(row['ingenannot_rna_aed'],'0.30000000' if r else '.')
                self.assertEqual(row['ingenannot_protein_aed'],'0.10000000' if p else '.')
                self.assertFalse(any('combined' in k for k in row))

    def test_changed_models_sources_and_ids_block_assessment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,c,i,b,a=self.fixtures(root)
            (root/'rna.gtf').write_text('# changed\n')
            with self.assertRaisesRegex(ValueError,'evidence changed'): check(gff,v,c,root/'rna.gtf',root/'protein.gff',root/'other.json')
            a.write_text(a.read_text().replace('exon\t101\t200','exon\t101\t199'))
            with self.assertRaisesRegex(ValueError,'changed model structure'): report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')

    def test_native_score_outside_range_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,c,i,b,a=self.fixtures(root)
            a.write_text(a.read_text().replace('aed_ev_tr=0.3000','aed_ev_tr=1.3000'))
            with self.assertRaisesRegex(ValueError,'Invalid native AED'): report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')

    def test_exon_envelope_normalization_preserves_blocks_on_both_strands(self):
        for strand in ('+','-'):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); gff,v,c,i,b,a=self.fixtures(root)
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
                data=report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')
                self.assertEqual(data['transcript_bounds_normalized_to_exons'],['t1'])
                a.write_text(a.read_text().replace('mRNA\t101','mRNA\t99'))
                with self.assertRaisesRegex(ValueError,'unexpected transcript bounds'):
                    report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')

    def test_cds_phase_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); gff,v,c,i,b,a=self.fixtures(root)
            lines=[]
            for line in a.read_text().splitlines():
                fields=line.split('\t')
                if len(fields)==9 and fields[2]=='CDS':
                    fields[7]='1' if fields[7]!='1' else '2'
                    line='\t'.join(fields)
                lines.append(line)
            a.write_text('\n'.join(lines)+'\n')
            with self.assertRaisesRegex(ValueError,'coding phases differ'):
                report(gff,a,i,root/'inputs.json',b,True,.25,root/'out')


if __name__=='__main__': unittest.main()
