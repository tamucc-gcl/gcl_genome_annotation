import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import subprocess

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from model_metrics import prepare
from model_assessment import report, load_table
from te_evidence import flag, repeat_bed
from te_search import search
from prepare_protein_evidence import checksum


def dump(path,value):path.write_text(json.dumps(value))


def fixtures(root):
    dna=list('A'*1000);dna[200:202]=list('GT');dna[298:300]=list('AG');dna[650:652]=list('CT');dna[748:750]=list('AC')
    genome=root/'genome.fa';genome.write_text('>chr1\n'+''.join(dna)+'\n')
    identity=root/'identity.json';dump(identity,dict(reference_identity='ref',contigs=[dict(id='chr1',length=1000,uppercase_sha256=hashlib.sha256(''.join(dna).encode()).hexdigest())]))
    gff=root/'A.gff3';lines=[]
    for gid,strand,bounds,tids,exons,cds in [('g1','+',(100,400),['t1','t2'],[(100,200),(300,400)],[(150,200),(300,350)]),
        ('g2','-',(600,800),['t3'],[(600,650),(750,800)],[(600,650),(750,800)])]:
        s,e=bounds;lines.append(f'chr1\tfixture\tgene\t{s+1}\t{e}\t.\t{strand}\t.\tID={gid};Note=old%3Bnote\n')
        for tid in tids:
            s,e=bounds
            lines.append(f'chr1\tfixture\tmRNA\t{s+1}\t{e}\t.\t{strand}\t.\tID={tid};Parent={gid}\n')
            for kind,blocks in [('exon',exons),('CDS',cds)]:
                for i,(s,e) in enumerate(blocks):lines.append(f'chr1\tfixture\t{kind}\t{s+1}\t{e}\t.\t{strand}\t{0 if kind=="CDS" else "."}\tID={tid}.{kind}{i};Parent={tid}\n')
    gff.write_text('##gff-version 3\n'+''.join(lines))
    proteins=root/'proteins.fa';proteins.write_text('>t1\nM'+'K'*32+'\n>t2\nMK*Y\n>t3\nM'+'P'*32+'*\n')
    cds=root/'cds.fa';cds.write_text(''.join(f'>{tid}\n'+ 'A'*100+'\n' for tid in ('t1','t2','t3')))
    state=root/'state.json';dump(state,dict(status='ready',stage='final_annotation',assembly_id='A',reference_identity='ref',
        gff3_sha256=checksum(gff),proteins_sha256=checksum(proteins),cds_sha256=checksum(cds)))
    metrics=root/'metrics';prepare('A',gff,proteins,cds,genome,identity,state,metrics)
    native=root/'search';native.mkdir();hits=native/'hits.tsv';hits.write_text('t1\trpt1\t90\t30\t33\t40\t1e-20\t100\t1\t30\n')
    dump(native/'search.json',dict(status='ready',assessed=True,reason=None,proteins_sha256=checksum(proteins),evalue=1e-10,hits_sha256=checksum(hits)))
    coverage=root/'coverage.tsv'
    with coverage.open('w') as h:
        for line in (metrics/'regions.bed').read_text().splitlines():
            r=line.split('\t');length=int(r[2])-int(r[1]);bp=length//2 if r[3].endswith('g1') else 0
            h.write(line+f'\t{int(bp>0)}\t{bp}\t{length}\t{bp/length}\n')
    return gff,proteins,cds,genome,identity,state,metrics,native,coverage


class ModelAssessmentContracts(unittest.TestCase):
    def test_reused_metrics_both_strands_isoforms_and_terminal_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root)
            data=json.loads((m/'metrics.json').read_text());tx={r['transcript_id']:r for r in data['metrics']}
            self.assertEqual(tx['t1']['donor_acceptor_summary'],'GTAG:1');self.assertEqual(tx['t3']['donor_acceptor_summary'],'GTAG:1')
            self.assertEqual(tx['t2']['has_internal_stop'],'Y');self.assertEqual(tx['t3']['has_internal_stop'],'N')
            self.assertEqual(tx['t3']['has_terminal_stop_symbol'],'Y')
            self.assertEqual(tx['t1']['cds_length_bp'],100)
            cds_regions=[r.split('\t') for r in (m/'regions.bed').read_text().splitlines() if 'cds|g1' in r]
            self.assertEqual(sum(int(r[2])-int(r[1]) for r in cds_regions),100)
            self.assertEqual(report('A',m/'metrics.json',root/'out')['models_removed'],0)

    def test_gene_flags_preserve_all_models_and_show_isoform_homology(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root);original=g.read_bytes()
            te=root/'te';state=flag('A',g,m/'metrics.json',coverage,n,.5,te)
            self.assertEqual(state['flagged_genes'],1)
            result=report('A',m/'metrics.json',root/'out',te_status=te/'te_status.json',te_scores=te/'te_gene_report.tsv')
            self.assertEqual(result['transcripts'],3);self.assertEqual(result['models_removed'],0);self.assertEqual(g.read_bytes(),original)
            rows=load_table(root/'out/transcript_assessment.tsv','transcript_id')
            self.assertEqual(rows['t1']['te_transcript_homology'],'True');self.assertEqual(rows['t2']['te_transcript_homology'],'False')
            self.assertEqual(rows['t2']['te_gene_flag'],'True');self.assertEqual(rows['t1']['rna_aed'],'.')
            flagged=(te/'A.te_flagged.gff3').read_text()
            self.assertIn('Note=old%3Bnote,putative_TE_derived',flagged)
            # Structure/IDs are preserved in a derived flag-only GFF.
            from evidence_aed import models
            self.assertEqual({k:v['cds'] for k,v in models(g,{'chr1':1000}).items()},
                             {k:v['cds'] for k,v in models(te/'A.te_flagged.gff3',{'chr1':1000}).items()})

    def test_no_hits_are_valid_but_failed_searches_are_not_suppressed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root)
            db=root/'db';db.mkdir();dump(db/'database.json',dict(status='ready',assessed=True,reason=None,library_sha256='lib'))
            with patch('te_search.subprocess.run',side_effect=subprocess.CalledProcessError(1,['diamond'])):
                with self.assertRaises(subprocess.CalledProcessError):search(p,db,1,1e-10,root/'bad')
            def success(command,check):Path(command[command.index('--out')+1]).write_text('')
            with patch('te_search.subprocess.run',side_effect=success):search(p,db,1,1e-10,root/'empty')
            status=flag('A',g,m/'metrics.json',coverage,root/'empty',.75,root/'te')
            self.assertEqual(status['flagged_genes'],0)

    def test_missing_provenance_wrong_ids_and_mutated_evidence_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root)
            p.write_text(p.read_text()+'# changed\n')
            with self.assertRaisesRegex(ValueError,'file changed'):prepare('A',g,p,c,genome,i,v,root/'bad')
            aed=root/'aed.tsv';aed.write_text('gene_id\ttranscript_id\trna_aed\tprotein_aed\ng1\tmissing\t.\t.\n')
            with self.assertRaisesRegex(ValueError,'validation status'):report('A',m/'metrics.json',root/'bad',aed_scores=aed)
            status=root/'aed_status.json';dump(status,dict(status='ready',model_sha256=checksum(g),reference_identity='ref'))
            with self.assertRaisesRegex(ValueError,'IDs differ'):report('A',m/'metrics.json',root/'bad',aed_status=status,aed_scores=aed)
            repeat=root/'repeat.gff';repeat.write_text('chr1\tRM\trepeat_region\t999\t1001\t.\t.\t.\tID=r1\n')
            with self.assertRaisesRegex(ValueError,'outside accepted'):repeat_bed(repeat,i,root/'bad.bed')

    def test_cluster_assessment_fixture_preserves_all_four_evidence_modes(self):
        from validate_assessment_cluster import fixtures as cluster_fixture
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cluster_fixture(root)
            state=json.loads((root/'accepted/structural_annotation/prediction_status.json').read_text())
            self.assertEqual({r['mode'] for r in state['assemblies']},{'rna_protein','rna','protein','ab_initio'})
            for row in state['assemblies']:
                aid=row['assembly_id'];folder=root/f'accepted/structural_annotation/final/{aid}'
                data=prepare(aid,folder/f'{aid}.gff3',folder/f'{aid}_proteins.fa',folder/f'{aid}_cds.fa',
                    root/f'accepted/repeats/{aid}.masked.fasta',root/f'accepted/filtered_assemblies/{aid}/{aid}.contigs.json',
                    folder/'annotation_status.json',root/f'metrics_{aid}')
                self.assertEqual(data['transcripts'],3)

    def test_empty_repeat_library_marks_homology_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root)
            db=root/'empty_db';db.mkdir();dump(db/'database.json',dict(status='ready',assessed=False,reason='Empty library',library_sha256='empty'))
            with patch('te_search.subprocess.run') as command:
                search(p,db,1,1e-10,root/'empty_search');command.assert_not_called()
            flag('A',g,m/'metrics.json',coverage,root/'empty_search',.75,root/'te')
            report('A',m/'metrics.json',root/'out',te_status=root/'te/te_status.json',te_scores=root/'te/te_gene_report.tsv')
            row=load_table(root/'out/transcript_assessment.tsv','transcript_id')['t1']
            self.assertEqual(row['te_transcript_homology'],'.');self.assertEqual(row['te_gene_flag'],'False')

    def test_accepted_evidence_join_and_reference_mismatch(self):
        from validate_assessment_cluster import fixtures as cluster_fixture
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cluster_fixture(root);accepted=root/'accepted'
            for aid in ('A','B','C','D'):
                final=accepted/f'structural_annotation/final/{aid}';aed=accepted/f'structural_annotation/aed/{aid}';j=accepted/f'qc/junctions/{aid}'
                m=root/f'metrics_{aid}'
                prepare(aid,final/f'{aid}.gff3',final/f'{aid}_proteins.fa',final/f'{aid}_cds.fa',accepted/f'repeats/{aid}.masked.fasta',
                    accepted/f'filtered_assemblies/{aid}/{aid}.contigs.json',final/'annotation_status.json',m)
                result=report(aid,m/'metrics.json',root/f'out_{aid}',aed_status=aed/'aed_status.json',aed_scores=aed/'aed_scores.tsv',
                    junction_status=j/'junction_status.json',junction_scores=j/'transcript_junction_support.tsv')
                self.assertEqual(result['transcripts'],3);self.assertEqual(result['models_removed'],0)
            path=accepted/'qc/junctions/A/junction_status.json';state=json.loads(path.read_text());state['reference_identity']='other';dump(path,state)
            with self.assertRaisesRegex(ValueError,'Evidence/model/reference mismatch'):
                report('A',root/'metrics_A/metrics.json',root/'bad',junction_status=path,
                    junction_scores=accepted/'qc/junctions/A/transcript_junction_support.tsv')

    def test_single_exon_metrics_are_not_zero_canonical_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,p,c,genome,i,v,m,n,coverage=fixtures(root)
            g.write_text(''.join(line for line in g.read_text().splitlines(keepends=True) if 'ID=t3.exon1;' not in line and 'ID=t3.CDS1;' not in line))
            state=json.loads(v.read_text());state['gff3_sha256']=checksum(g);dump(v,state)
            data=prepare('A',g,p,c,genome,i,v,root/'single')
            row=next(r for r in data['metrics'] if r['transcript_id']=='t3')
            self.assertEqual(row['n_introns'],0);self.assertIsNone(row['pct_canonical_splice']);self.assertIsNone(row['min_intron_len'])


if __name__=='__main__':unittest.main()
