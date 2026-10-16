from pathlib import Path
import csv
import json
import tempfile
import unittest

from review_ingenannot_evidence import compare_rna, selected_proteins, review
from prepare_protein_evidence import checksum


class CoordinateReviewContracts(unittest.TestCase):
    def test_review_reads_completed_outputs_and_writes_matching_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); s=root/'structural_annotation'
            gff=s/'braker/A/A.gff3'; gff.parent.mkdir(parents=True)
            gff.write_text(''.join(f'chr1\tBRAKER\t{k}\t101\t200\t.\t+\t{phase}\t{atts}\n'
                for k,phase,atts in [('gene','.','ID=g1'),('mRNA','.','ID=t1;Parent=g1'),
                    ('exon','.','Parent=t1'),('CDS','0','Parent=t1')]))
            rna=gff.parent/'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
            rna.parent.mkdir(parents=True)
            rna.write_text(''.join(f'chr1\tStringTie\t{k}\t101\t200\t.\t+\t.\tgene_id "r"; transcript_id "r1";\n' for k in ('transcript','exon')))
            identity=root/'filtered_assemblies/A/A.contigs.json'; identity.parent.mkdir(parents=True)
            identity.write_text(json.dumps(dict(reference_identity='ref',contigs=[dict(id='chr1',length=1000)])))
            bench=s/'ingenannot_benchmark/A'; bench.mkdir(parents=True)
            (bench/'benchmark_status.json').write_text(json.dumps(dict(status='ready',model_sha256=checksum(gff),reference_identity='ref',
                profile=dict(rna_cds_only=True,rna_overflow_penalty=.25),evidence=dict(rna=dict(source_sha256=checksum(rna))))))
            evidence=s/'evidence_aed/A'; (evidence/'protein_alignment').mkdir(parents=True)
            (evidence/'protein_alignment/A.miniprot.gff').write_text('chr1\tminiprot\tmRNA\t101\t200\t.\t+\t.\tID=MP1;Target=p 1 33\nchr1\tminiprot\tCDS\t101\t200\t.\t+\t0\tParent=MP1\n')
            row=dict(gene_id='g1',transcript_id='t1',current_rna_aed='0',ingenannot_rna_aed='0',rna_delta='0',
                rna_penalty='no',rna_evidence_id='r1',current_protein_aed='0',ingenannot_protein_aed='0',protein_evidence_id='MP1')
            with (bench/'benchmark_scores.tsv').open('w',newline='') as handle:
                w=csv.DictWriter(handle,fieldnames=list(row),delimiter='\t'); w.writeheader(); w.writerow(row)
            (evidence/'aed_scores.tsv').write_text('transcript_id\trna_aed\tprotein_aed\trna_best_alignment_id\nt1\t0\t0\tr1\n')
            result=review(root,'A',root/'out')
            self.assertEqual(result['status'],'passed')
            self.assertEqual(result['cases'][0]['category'],'strong_agreement_control')
            self.assertTrue((root/'out/review.json').is_file())

    def test_utr_overhangs_do_not_trigger_penalty_but_coding_conflicts_do(self):
        model=dict(cds=[(100,150),(200,250)],exons=[(90,150),(200,270)])
        rna=dict(identifier='r1',seq='chr1',strand='+',blocks=[(80,150),(200,280)])
        data=compare_rna(model,rna,.25)
        self.assertFalse(data['penalty_applied'])
        self.assertEqual(data['matching_intron_boundaries'],1)
        rna['blocks']=[(80,151),(200,280)]
        data=compare_rna(model,rna,.25)
        self.assertTrue(data['penalty_applied'])
        self.assertEqual(data['rna_extra_exonic_bases_inside_cds_span'],1)
        self.assertEqual(data['matching_intron_boundaries'],0)
        self.assertAlmostEqual(data['penalized_aed']-data['raw_cds_aed'],.25)
        rna['blocks']=[(80,149),(200,280)]
        self.assertEqual(compare_rna(model,rna,.25)['cds_bases_outside_rna'],1)

    def test_selected_protein_scan_keeps_all_blocks_and_exact_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'proteins.gff'
            p.write_text('##gff-version 3\n'+''.join(
                f'chr1\tminiprot\t{kind}\t{s}\t{e}\t.\t+\t.\t{attrs}\n'
                for kind,s,e,attrs in [('mRNA',1,30,'ID=MP10;Target=p 1 10'),
                    ('CDS',1,30,'Parent=MP10'),('mRNA',101,250,'ID=MP1;Target=q 1 30'),
                    ('CDS',101,150,'Parent=MP1'),('CDS',201,250,'Parent=MP1'),
                    ('mRNA',300,350,'ID=MP2;Target=x 1 10')]))
            hits=selected_proteins(p,{'MP1'})
            self.assertEqual(set(hits),{'MP1'})
            self.assertEqual(hits['MP1']['blocks'],[(100,150),(200,250)])
            with self.assertRaisesRegex(ValueError,'missing'): selected_proteins(p,{'MP999'})


if __name__=='__main__': unittest.main()
