import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from model_selection import select, settings, table, write_table
from model_assessment import report
from model_metrics import prepare, sequences
from annotation_models import models
from prepare_protein_evidence import checksum
from test_model_assessment import fixtures


def fixture(root,single=False,shared=False):
    g,p,c,genome,i,v,m,n,coverage=fixtures(root)
    if single:
        g.write_text(''.join(line for line in g.read_text().splitlines(keepends=True) if not any(x in line for x in ('ID=t3.exon1;','ID=t3.CDS1;'))))
    if shared:
        lines=[]
        for line in g.read_text().splitlines(keepends=True):
            if 'Parent=t2' in line:continue
            if 'Parent=t1' in line:line=line.replace('ID=t1.','ID=shared.').replace('Parent=t1','Parent=t1,t2')
            lines.append(line)
        g.write_text(''.join(lines))
    state=json.loads(v.read_text());state['gff3_sha256']=checksum(g);v.write_text(json.dumps(state))
    prepare('A',g,p,c,genome,i,v,m)
    folder=root/'assessment';report('A',m/'metrics.json',folder)
    return g,p,c,i,v,folder


def change_assessment(folder,changes,**flags):
    path=folder/'transcript_assessment.tsv';fields,rows=table(path,'transcript_id')
    for row in rows:row.update(changes.get(row['transcript_id'],{}))
    write_table(path,fields,rows)
    status=folder/'assessment_status.json';state=json.loads(status.read_text());state.update(flags)
    state['outputs']['transcript_assessment.tsv']=checksum(path);status.write_text(json.dumps(state))


def run(root,inputs,options=None):
    g,p,c,i,v,folder=inputs
    return select('A',g,p,c,v,i,folder/'transcript_assessment.tsv',folder/'assessment_status.json',root/'selection',options)


class ModelSelectionContracts(unittest.TestCase):
    def test_original_defaults_and_exact_sequence_partition(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,single=True);state=run(root,data)
            self.assertEqual((state['retained_transcripts'],state['excluded_transcripts'],state['excluded_genes']),(2,1,1))
            out=root/'selection'
            self.assertEqual(set(sequences(out/'retained/A_proteins.fa')),{'t1','t2'})
            self.assertEqual(set(sequences(out/'excluded/A_cds.fa')),{'t3'})
            self.assertEqual(state['removal_reason_counts'],{'single_exon_short_cds':1})
            self.assertEqual(sequences(data[1])['t1'],sequences(out/'retained/A_proteins.fa')['t1'])
            self.assertTrue(all(checksum(out/name)==sha for name,sha in state['outputs'].items()))

    def test_shared_features_split_and_good_isoform_rescues_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,shared=True);folder=data[-1]
            change_assessment(folder,{'t1':{'protein_aed':'.2'},'t2':{'protein_aed':'.8'},'t3':{'protein_aed':'.1'}},aed_enabled=True)
            state=run(root,data,dict(max_aed=.5,aed_source='protein'))
            self.assertEqual(state['retained_genes'],2);self.assertEqual(state['excluded_transcripts'],1)
            for label,tid in [('retained','t1'),('excluded','t2')]:
                text=(root/f'selection/{label}/A.gff3').read_text()
                self.assertIn('ID=shared.exon0;Parent='+tid,text);self.assertNotIn('Parent=t1,t2',text)
            self.assertEqual(set(models(root/'selection/excluded/A.gff3',{'chr1':1000})),{'t2'})

    def test_boundary_and_flag_only_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,single=True)
            state=run(root,data,dict(min_cds_length=50,min_gene_length=200,min_intron_length=100))
            self.assertEqual(state['excluded_transcripts'],0) # Equality retains at every structural boundary.
            state=run(root,data,dict(structural_action='flag'))
            self.assertEqual(state['excluded_transcripts'],0)
            rows=table(root/'selection/transcript_decisions.tsv','transcript_id')[1]
            self.assertIn('single_exon_short_cds:flag',next(r for r in rows if r['transcript_id']=='t3')['rule_hits'])
            self.assertEqual(run(root,data,dict(structural_action='off'))['excluded_transcripts'],0)

    def test_no_rna_and_single_exon_junctions_are_not_zero_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,single=True)
            change_assessment(data[-1],{'t1':{'junction_coding_pass_fraction':'.9'},
                't2':{'junction_coding_pass_fraction':'.','junction_coding_junctions':'1'},
                't3':{'junction_coding_pass_fraction':'.','junction_coding_junctions':'0'}},junction_qc_enabled=True)
            state=run(root,data,dict(structural_action='off',min_junction_fraction=.9))
            self.assertEqual(state['excluded_transcripts'],0)
            rows={r['transcript_id']:r for r in table(root/'selection/transcript_decisions.tsv','transcript_id')[1]}
            self.assertIn('coding_junction_support_unavailable',rows['t2']['unavailable_or_not_applicable'])
            self.assertIn('coding_junctions_not_applicable',rows['t3']['unavailable_or_not_applicable'])

    def test_all_triggered_reasons_te_gene_level_and_empty_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root,single=True)
            change_assessment(data[-1],{t:dict(te_gene_flag='True',protein_aed='.8') for t in ['t1','t2','t3']},te_enabled=True,aed_enabled=True)
            state=run(root,data,dict(te_action='remove',max_aed=.2,aed_source='protein'))
            self.assertTrue(state['empty_retained_set']);self.assertEqual(state['excluded_genes'],2)
            rows=table(root/'selection/transcript_decisions.tsv','transcript_id')[1]
            self.assertEqual(set(next(r for r in rows if r['transcript_id']=='t3')['removal_reasons'].split(';')),
                {'single_exon_short_cds','high_protein_aed','gene_repeat_association'})
            self.assertEqual(sequences(root/'selection/retained/A_cds.fa'),{})

    def test_unavailable_aed_skipped_disabled_stage_blocks_requested_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root)
            with self.assertRaisesRegex(ValueError,'requires enabled AED'):run(root,data,dict(max_aed=.5,aed_source='rna'))
            change_assessment(data[-1],{},aed_enabled=True)
            state=run(root,data,dict(max_aed=.5,aed_source='rna'))
            self.assertEqual(state['excluded_transcripts'],0)
            self.assertTrue(all(r['unavailable_or_not_applicable'].startswith('rna_aed_unavailable') for r in table(root/'selection/transcript_decisions.tsv','transcript_id')[1]))

    def test_mutated_sources_and_wrong_parents_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root)
            change_assessment(data[-1],{'t1':{'gene_id':'wrong'}})
            with self.assertRaisesRegex(ValueError,'gene parent'):run(root,data)
            data[1].write_text(data[1].read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'sequences changed'):run(root,data)

    def test_invalid_parameters_and_scientific_notation(self):
        self.assertEqual(settings(dict(max_aed='1e-2',aed_source='rna'))['max_aed'],.01)
        for opts in [dict(max_aed=.5),dict(min_cds_length=-1),dict(min_gene_length=1.5),dict(te_action='yes'),dict(max_aed='NaN'),dict(min_junction_fraction=2),dict(unknown=1)]:
            with self.assertRaises(ValueError):settings(opts)

    def test_short_introns_inferred_without_intron_features_and_gene_span_rule(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);data=fixture(root)
            self.assertNotIn('\tintron\t',data[0].read_text())
            state=run(root,data,dict(min_intron_length=101))
            self.assertEqual(state['excluded_transcripts'],3)
            self.assertEqual(state['removal_reason_counts'],{'all_introns_short':3})
            state=run(root,data,dict(min_gene_length=250))
            self.assertEqual(state['excluded_transcripts'],1)
            self.assertEqual(state['removal_reason_counts'],{'short_gene_span':1})


if __name__=='__main__':unittest.main()
