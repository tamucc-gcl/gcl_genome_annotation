import argparse
import base64
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from functional_inputs import TOOLS, prepare, read_table
from merge_functional_annotations import merge
from annotate_gff3 import annotate
from functional_report import report, geometry
from prepare_protein_evidence import checksum


def fixture(root, aid='A'):
    folder=root/aid;folder.mkdir(parents=True)
    gff=folder/f'{aid}.gff3';proteins=folder/f'{aid}_proteins.fa';assessment=folder/'transcript_assessment.tsv'
    rows=[]
    for i,strand in [(1,'+'),(2,'-'),(3,'+')]:
        tid=f'{aid}_g{i}.t1';gid=f'{aid}_g{i}';start=i*100
        for kind,attrs,phase in [('gene',f'ID={gid}','.'),('mRNA',f'ID={tid};Parent={gid}','.'),
                                 ('exon',f'ID={tid}.exon1;Parent={tid}','.'),('CDS',f'ID={tid}.CDS1;Parent={tid}','0')]:
            rows.append(f'chr1\toriginal\t{kind}\t{start}\t{start+59}\t.\t{strand}\t{phase}\t{attrs}\n')
    gff.write_text('##gff-version 3\n'+''.join(rows))
    proteins.write_text(f'>{aid}_g1.t1\nMALWMRLLPLLALLALWGPGPG*\n>{aid}_g2.t1\nALWMRLLPLLALLALWGPGPG\n>{aid}_g3.t1\nMA*AA\n')
    assessment.write_text('transcript_id\tgene_id\tstructure_flag\trna_aed\tprotein_aed\tjunction_status\tte_gene_flag\n'+
        ''.join(f'{aid}_g{i}.t1\t{aid}_g{i}\t'+('FAIL' if i==3 else 'PASS')+'\t0.2\t0.1\tready\tfalse\n' for i in (1,2,3)))
    validation=folder/'annotation_status.json'
    validation.write_text(json.dumps(dict(status='ready',stage='identified_annotation',assembly_id=aid,coding_transcripts=3,
        gff3_sha256=checksum(gff),proteins_sha256=checksum(proteins),outputs={'transcript_assessment.tsv':checksum(assessment)})))
    return folder


def arguments(root, aid='A', enabled=True):
    folder=fixture(root,aid);prep=root/'prepared'
    prepare(aid,folder/f'{aid}_proteins.fa',folder/f'{aid}.gff3',folder/'transcript_assessment.tsv',folder/'annotation_status.json',prep)
    paths={}
    for name in ('interproscan_tsv','eggnog_tsv','diamond_tsv','kofamscan_tsv','deepsig_gff3'):
        path=root/name;path.write_text('');paths[name]=str(path)
    settings={f'run_{tool}':enabled for tool in TOOLS}
    return argparse.Namespace(assembly_id=aid,assessment=str(folder/'transcript_assessment.tsv'),eligibility=str(prep/'eligibility.tsv'),
        settings=base64.b64encode(json.dumps(settings).encode()).decode(),omark=None,outdir=str(root/'merged'),**paths),folder


def decorate(args,folder,root):
    options=argparse.Namespace(assembly_id=args.assembly_id,braker_gff3=str(folder/f'{args.assembly_id}.gff3'),
        functional_tsv=str(root/f'merged/{args.assembly_id}.functional_annotations.tsv'),
        **{k:os.devnull for k in ('trnascan_gff3','infernal_gff3','telomere_gff3','centromere_gff3','repeat_gff')})
    previous=Path.cwd()
    try:
        os.chdir(root);annotate(options)
    finally:os.chdir(previous)


class FunctionalContracts(unittest.TestCase):
    def test_eligibility_preserves_partial_and_does_not_delete_internal_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root)
            prepared=(root/'prepared/eligible.fa').read_text()
            self.assertIn('A_g2.t1',prepared);self.assertNotIn('A_g3.t1',prepared);self.assertNotIn('*',prepared)
            self.assertIn('MA*AA',(folder/'A_proteins.fa').read_text())

    def test_empty_outputs_are_no_hit_and_disabled_are_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root)
            state=merge(args);self.assertEqual(state['tools']['diamond'],{'no_hit':2,'ineligible':1})
            decorate(args,folder,root)
            self.assertIn('product=hypothetical protein',(root/'A.annotated.gff3').read_text())
            state=report('A',folder/'A.gff3',root/'A.annotated.gff3',root/'merged/A.functional_annotations.tsv',
                root/'merged/functional_merge_status.json',args.eligibility,root/'report')
            self.assertEqual(state['hypothetical_transcripts'],3)
            args.settings=base64.b64encode(json.dumps({f'run_{t}':False for t in TOOLS}).encode()).decode()
            state=merge(args);self.assertEqual(state['tools']['diamond'],{'disabled':3})

    def test_native_formats_product_hierarchy_and_separate_evidence_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root);tid='A_g1.t1'
            Path(args.interproscan_tsv).write_text('\t'.join([tid,'md5','22','Pfam','PF0001','domain','1','20','1e-9','T','date','IPR0001','desc','GO:0000001'])+'\n')
            egg=[tid,'seed','1e-20','100','OG1@1','1','A','egg description','eggname','GO:0000002','-','ko:K00001','map00010','-','-','GH1','-','-','-','-','PF0001']
            Path(args.eggnog_tsv).write_text('#header\n'+'\t'.join(egg)+'\n')
            Path(args.diamond_tsv).write_text('qseqid\tsseqid\tstitle\tpident\tevalue\n'+tid+'\tsp|P12345|TEST\tsp|P12345|TEST Good; protein, name OS=Species\t99\t1e-30\n')
            Path(args.kofamscan_tsv).write_text('# marker\tprotein\tKO\tthreshold\tscore\tevalue\tdefinition\n'+
                f'*\t{tid}\tK00001\t10\t100\t1e-20\tKO name\n\tA_g2.t1\tK00002\t100\t1\t0.1\tNot significant\n')
            Path(args.deepsig_gff3).write_text(f'{tid}\tDeepSig\tSignal peptide\t1\t10\t0.9\t.\t.\tevidence=x\nA_g2.t1\tDeepSig\tChain\t1\t21\t.\t.\t.\tevidence=x\n')
            state=merge(args);decorate(args,folder,root)
            _,rows=read_table(root/'merged/A.functional_annotations.tsv','protein_id')
            self.assertEqual(rows['A_g2.t1']['kofamscan_status'],'no_hit');self.assertEqual(rows['A_g2.t1']['deepsig_status'],'no_hit')
            self.assertEqual(rows[tid]['rna_aed'],'0.2');self.assertEqual(rows[tid]['protein_aed'],'0.1')
            self.assertEqual(rows[tid]['structure_flag'],'PASS');self.assertEqual(rows[tid]['junction_status'],'ready')
            text=(root/'A.annotated.gff3').read_text();self.assertIn('Good%3B protein%2C name',text)
            self.assertNotIn('_AED=',text);self.assertIn('rna_aed=0.2;protein_aed=0.1',text)
            self.assertEqual(geometry(folder/'A.gff3'),geometry(root/'A.annotated.gff3'))

    def test_unknown_and_ineligible_hit_ids_are_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root)
            for tid in ('foreign','A_g3.t1'):
                Path(args.diamond_tsv).write_text('qseqid\tsseqid\n'+tid+'\tx\n')
                with self.assertRaises(ValueError):merge(args)

    def test_source_change_and_coordinate_change_are_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root);merge(args);decorate(args,folder,root)
            target=root/'A.annotated.gff3';target.write_text(target.read_text().replace('\t100\t159','\t101\t159'))
            with self.assertRaisesRegex(ValueError,'geometry'):
                report('A',folder/'A.gff3',target,root/'merged/A.functional_annotations.tsv',root/'merged/functional_merge_status.json',args.eligibility,root/'report')
            proteins=folder/'A_proteins.fa';proteins.write_text(proteins.read_text()+'A')
            with self.assertRaisesRegex(ValueError,'changed'):
                prepare('A',proteins,folder/'A.gff3',folder/'transcript_assessment.tsv',folder/'annotation_status.json',root/'other')

    def test_omark_classifies_representatives_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root);ump=root/'A.ump';ump.write_text('>Consistent\nA_g1.t1\n');args.omark=str(ump)
            merge(args);_,rows=read_table(root/'merged/A.functional_annotations.tsv','protein_id')
            self.assertEqual(rows['A_g1.t1']['omark_status'],'classified')
            self.assertEqual(rows['A_g2.t1']['omark_status'],'not_classified')

    def test_columns_remain_stable_and_below_threshold_ko_does_not_name_a_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root)
            Path(args.kofamscan_tsv).write_text('\tA_g1.t1\tK00001\t100\t10\t0.1\tBelow threshold name\n')
            merge(args);columns,_=read_table(root/'merged/A.functional_annotations.tsv','protein_id')
            self.assertIn('diamond_hit',columns);self.assertIn('kofam_KO',columns)
            decorate(args,folder,root)
            self.assertNotIn('Below threshold name',(root/'A.annotated.gff3').read_text())
            self.assertIn('product=hypothetical protein',(root/'A.annotated.gff3').read_text())

    def test_best_significant_ko_is_deterministic_and_malformed_output_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args,folder=arguments(root)
            Path(args.kofamscan_tsv).write_text('*\tA_g1.t1\tK00001\t10\t100\t1e-20\tBest KO\n*\tA_g1.t1\tK00002\t10\t20\t1e-5\tOther KO\n')
            merge(args);_,rows=read_table(root/'merged/A.functional_annotations.tsv','protein_id')
            self.assertEqual(rows['A_g1.t1']['kofam_KO'],'K00001')
            Path(args.interproscan_tsv).write_text('A_g1.t1\ttruncated\n')
            with self.assertRaisesRegex(ValueError,'Malformed'):merge(args)


if __name__=='__main__':unittest.main()
