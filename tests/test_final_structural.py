import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from final_identifiers import rename
from structural_report import report
from model_selection import table, select
from model_metrics import sequences
from evidence_aed import models
from prepare_protein_evidence import checksum
from annotation_baseline import prepare
from identified_qc import prepare_selected
from test_model_selection import fixture,change_assessment,run


def accepted(root):
    data=fixture(root,single=True,shared=True)
    g,p,c,identity,validation,assessment=data
    # Real content-derived identity, as required by native QC representative preparation.
    ref=json.loads(identity.read_text());ref['reference_identity']=hashlib.sha256(json.dumps(ref['contigs'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    identity.write_text(json.dumps(ref))
    for path in (validation,assessment/'assessment_status.json'):
        state=json.loads(path.read_text());state['reference_identity']=ref['reference_identity'];path.write_text(json.dumps(state))
    # Model mode does not affect selection/final IDs.
    run(root,data)
    final=root/'identified';rename('A',root/'selection',identity,final)
    return data,final


def qc_fixture(root,final):
    meta=json.loads((final/'annotation_status.json').read_text());qc=root/'qc';qc.mkdir(exist_ok=True)
    settings=dict(run_busco=True,run_omark=True,busco_lineage='test_odb10',taxid='1')
    prepare('A',settings,meta,final/'A.gff3',final/'A_proteins.fa',qc)
    (qc/'A_annotation_stats').write_text('Number of gene 1\nNumber of mRNA 2\n')
    (qc/'A_busco.txt').write_text('C:90.0%[S:85.0%,D:5.0%],F:5.0%,M:5.0%,n:100\n')
    (qc/'A.summary.txt').write_text('S:80.0%,D:10.0%[U:5.0%,E:5.0%],M:10.0%\nA:90.0%[P:1.0%,F:2.0%],I:4.0%,C:1.0%,U:5.0%\n')
    ids=list(sequences(final/'A_proteins.fa'));(qc/'A.ump').write_text('>Consistent\n'+ids[0]+'\n')
    return [qc/n for n in ['selection.json','representatives.tsv','A_annotation_stats','A_busco.txt','A.summary.txt','A.ump']]


class FinalStructuralContracts(unittest.TestCase):
    def test_representative_ties_follow_original_ids_not_renumbering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,final=accepted(root)
            p=final/'A_proteins.fa';p.write_text('>A_g1.t1\nMKKK\n>A_g1.t2\nMPPP\n')
            state=json.loads((final/'annotation_status.json').read_text());state['proteins_sha256']=checksum(p)
            state['representative_tie_ids']={'A_g1.t1':'z_original','A_g1.t2':'a_original'}
            result=prepare_selected('A',dict(run_busco=False,run_omark=False),state,final/'A.gff3',p,root/'qc')
            self.assertEqual(set(sequences(root/'qc/A.fa')),{'A_g1.t2'})
            self.assertEqual(result['model_set'],'selected')

    def test_complete_maps_sequences_phase_and_excluded_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,final=accepted(root)
            state=json.loads((final/'annotation_status.json').read_text())
            self.assertEqual((state['genes'],state['coding_transcripts']),(1,2))
            tx=models(final/'A.gff3',{'chr1':1000});self.assertEqual(set(tx),{'A_g1.t1','A_g1.t2'})
            self.assertEqual(sequences(final/'A_cds.fa'),{'A_g1.t1':'A'*100,'A_g1.t2':'A'*100})
            fields,rows=table(final/'transcript_id_map.tsv','transcript_id');mapping={r['transcript_id']:r['final_transcript_id'] for r in rows}
            self.assertEqual(mapping,dict(t1='A_g1.t1',t2='A_g1.t2',t3='.'))
            _,features=table_map(final/'feature_id_map.tsv')
            self.assertTrue(any(r['old_parent']=='t1,t2' and r['new_parent']=='A_g1.t1,A_g1.t2' for r in features))
            self.assertTrue(any(r['partition']=='excluded' and r['new_id']=='.' for r in features))
            _,rows=table(final/'transcript_assessment.tsv','transcript_id')
            self.assertEqual({r['original_transcript_id'] for r in rows},{'t1','t2'})

    def test_reordered_features_keep_id_assignments_and_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,single=True,shared=True)
            # Input gene order differs from reference/coordinate order, while children precede parents.
            g=data[0];g.write_text('##gff-version 3\n'+''.join(reversed(g.read_text().splitlines(keepends=True)[1:])))
            state=json.loads(data[4].read_text());state['gff3_sha256']=checksum(g);data[4].write_text(json.dumps(state))
            from model_metrics import prepare as metrics
            metrics('A',data[0],data[1],data[2],root/'genome.fa',data[3],data[4],root/'metrics')
            from model_assessment import report as assessment
            assessment('A',root/'metrics/metrics.json',data[-1]);run(root,data)
            final=root/'out';rename('A',root/'selection',data[3],final,'fish_A')
            self.assertEqual(set(sequences(final/'A_proteins.fa')),{'fish_A_g1.t1','fish_A_g1.t2'})

    def test_multipart_and_anonymous_features_have_complete_maps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,shared=True)
            g=data[0];text=g.read_text().replace('ID=shared.CDS1','ID=shared.CDS0').replace('ID=t3.exon1;','')
            g.write_text(text);state=json.loads(data[4].read_text());state['gff3_sha256']=checksum(g);data[4].write_text(json.dumps(state))
            from model_metrics import prepare as metrics
            from model_assessment import report as assessment
            metrics('A',g,data[1],data[2],root/'genome.fa',data[3],data[4],root/'metrics');assessment('A',root/'metrics/metrics.json',data[-1]);run(root,data)
            final=root/'out';rename('A',root/'selection',data[3],final)
            _,rows=table_map(final/'feature_id_map.tsv')
            repeated=[r for r in rows if r['old_id']=='shared.CDS0' and r['partition']=='retained']
            self.assertEqual(len(repeated),2);self.assertEqual(len({r['new_id'] for r in repeated}),1)
            self.assertTrue(any(r['old_id']=='.' and r['new_id']!='.' for r in rows))

    def test_mutated_input_invalid_prefix_and_empty_selection_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,final=accepted(root)
            with self.assertRaisesRegex(ValueError,'prefix'):rename('A',root/'selection',data[3],root/'bad','bad prefix')
            p=root/'selection/retained/A_cds.fa';p.write_text(p.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'changed'):rename('A',root/'selection',data[3],root/'bad')
            run(root,data,dict(min_gene_length=10000))
            with self.assertRaisesRegex(ValueError,'No retained'):rename('A',root/'selection',data[3],root/'bad')

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'Report figures verified in the plotting environment/native cluster task')
    def test_report_native_summaries_final_ids_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,final=accepted(root);qc=qc_fixture(root,final)
            state=report('A',final,root/'selection',data[-1]/'transcript_assessment.tsv',data[-1]/'assessment_status.json',[],qc,root/'report')
            self.assertEqual(state['model_sets'][1]['genes'],1);self.assertFalse(state['confidence_list_created'])
            _,metrics=table(root/'report/completeness_metrics.tsv','metric')
            self.assertEqual(next(r['value'] for r in metrics if r['metric']=='busco_complete'),'90.0')
            self.assertEqual(len(list((root/'report/plots').glob('*.svg'))),10)
            self.assertIn('Functional annotation remains pending',(root/'report/report.html').read_text())
            _,rows=table(root/'report/structural_qc_transcripts.tsv','transcript_id')
            self.assertTrue(all(r['gene_id']=='A_g1' for r in rows));self.assertEqual(rows[0]['omark_representative_classification'],'Consistent')

    def test_wrong_native_qc_model_and_unrecognized_summary_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data,final=accepted(root);qc=qc_fixture(root,final)
            meta=qc[0];value=json.loads(meta.read_text());value['gff3_sha256']='wrong';meta.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError,'different annotation'):report('A',final,root/'selection',data[-1]/'transcript_assessment.tsv',data[-1]/'assessment_status.json',[],qc,root/'bad')
            value['gff3_sha256']=checksum(final/'A.gff3');meta.write_text(json.dumps(value));qc[3].write_text('unexpected tool output')
            with self.assertRaisesRegex(ValueError,'Unrecognized BUSCO'):report('A',final,root/'selection',data[-1]/'transcript_assessment.tsv',data[-1]/'assessment_status.json',[],qc,root/'bad')


def table_map(path):
    import csv
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t');return reader.fieldnames,list(reader)


if __name__=='__main__':unittest.main()
