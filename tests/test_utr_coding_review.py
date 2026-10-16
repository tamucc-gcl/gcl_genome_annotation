import csv
import json
from pathlib import Path
import tempfile
import unittest

from prepare_protein_evidence import checksum
from review_utr_coding_evidence import review


class UtrCodingReview(unittest.TestCase):
    def test_existing_protein_chains_distinguish_neighboring_coding_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'structural_annotation'
            gff=base/'braker/A/A.gff3';gff.parent.mkdir(parents=True)
            rows=[]
            for n,start,end in [(1,1,30),(2,41,70)]:
                for kind,phase,atts in [('gene','.',f'ID=g{n}'),('mRNA','.',f'ID=t{n};Parent=g{n}'),
                    ('exon','.',f'ID=e{n};Parent=t{n}'),('CDS','0',f'ID=c{n};Parent=t{n}')]:
                    rows.append(f'chr1\tx\t{kind}\t{start}\t{end}\t.\t+\t{phase}\t{atts}\n')
            gff.write_text(''.join(rows))
            (gff.parent/'A_proteins.fa').write_text('>t1\nMMMMMMMMMM\n>t2\nAAAAAAAAAA\n')
            identity=root/'filtered_assemblies/A/A.contigs.json';identity.parent.mkdir(parents=True)
            identity.write_text(json.dumps(dict(reference_identity='ref',contigs=[dict(id='chr1',length=100)])))
            bench=base/'ingenannot_benchmark/A';bench.mkdir(parents=True)
            (bench/'benchmark_status.json').write_text(json.dumps(dict(status='ready',model_sha256=checksum(gff),reference_identity='ref')))
            with (bench/'benchmark_scores.tsv').open('w',newline='') as h:
                w=csv.DictWriter(h,fieldnames=['transcript_id','protein_evidence_id','ingenannot_protein_aed','ingenannot_rna_aed'],delimiter='\t')
                w.writeheader();w.writerows([dict(transcript_id=f't{i}',protein_evidence_id=f'MP{i}',
                    ingenannot_protein_aed='0',ingenannot_rna_aed='0.1') for i in (1,2)])
            alignment=base/'evidence_aed/A/protein_alignment/A.miniprot.gff';alignment.parent.mkdir(parents=True)
            alignment.write_text(''.join(f'chr1\tminiprot\tmRNA\t{s}\t{e}\t.\t+\t.\tID=MP{i};Target=protein{i} 1 10\n'
                f'chr1\tminiprot\tCDS\t{s}\t{e}\t.\t+\t0\tParent=MP{i}\n' for i,s,e in [(1,1,30),(2,41,70)]))
            selected=root/'review.json';selected.write_text(json.dumps(dict(selected_cases=[dict(transcript_id='t1',
                original=dict(seq='chr1',strand='+',cds=[[0,30]]),other_gene_coding_spans_overlapped=[dict(
                    gene_id='g2',strand='+',extension_overlap_cds_bp=30)])])))
            result=review(root,'A',selected,root/'out')
            self.assertEqual(len(result['models']),2)
            self.assertEqual(result['models'][0]['selected_protein_chain']['attributes']['Target'],'protein1 1 10')
            self.assertEqual(result['models'][0]['selected_chain_overlap_other_coding_models']['t2']['model_cds_fraction'],0)
            self.assertEqual(result['models'][1]['recomputed_protein_aed'],0)
