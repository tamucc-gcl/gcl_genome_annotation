import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from test_final_structural import accepted
from package_annotation import package, gff_records, extract_spliced, trna_models
from model_metrics import sequences
from reference_repeats import sha


def release_fixture(root):
    data, identified = accepted(root)
    identity = data[3]; ref = json.loads(identity.read_text())['reference_identity']
    def status(name, value):
        p = root/name; p.write_text(json.dumps(value)); return str(p)
    repeat = root/'repeat.gff3'
    repeat.write_text('##gff-version 3\nchr1\tRepeatMasker\trepeat_region\t1\t20\t.\t+\t.\tID=repeat1;Target=external_repeat 1 20 +\n')
    tracks = []; states = {}
    for t in ('trnascan','infernal','telomere','centromere'):
        if t == 'centromere': states[t] = dict(status='disabled'); continue
        p = root/f'A.{t}.gff3'; tracks.append(str(p))
        text = {'trnascan':'chr1\ttRNAscan\texon\t10\t30\t.\t+\t.\tID=A.trna.exon;Parent=A.trna\nchr1\ttRNAscan\ttRNA\t10\t30\t.\t+\t.\tID=A.trna\n'
                'chr1\ttRNAscan\ttRNA\t40\t70\t.\t-\t.\tID=A.trna.minus\n'
                'chr1\ttRNAscan\texon\t40\t45\t.\t-\t.\tID=A.trna.minus.e1;Parent=A.trna.minus\n'
                'chr1\ttRNAscan\texon\t60\t70\t.\t-\t.\tID=A.trna.minus.e2;Parent=A.trna.minus\n',
                'infernal':'', 'telomere':'chr1\ttidk\ttelomeric_repeat\t1\t10\t.\t-\t.\tID=A.telomere1\n'}[t]
        p.write_text('##gff-version 3\n'+text)
        states[t] = dict(status='ready',reference_identity=ref,features=len(text.splitlines()),gff3_sha256=sha(p))
    track_status = status('tracks.json',dict(status='ready',stage='genome_tracks',assembly_id='A',reference_identity=ref,tracks=states))
    repeat_status = status('repeats.json',dict(status='ready',stage='reference_repeats',assembly_id='A',reference_identity=ref,masked_sha256=sha(root/'genome.fa')))
    masked = root/'masked.fa'; masked.write_bytes((root/'genome.fa').read_bytes())
    coding = root/'decorated.gff3'; coding.write_text((identified/'A.gff3').read_text().replace('ID=A_g1;', 'ID=A_g1;product=test;'))
    table = root/'functional.tsv'; table.write_text('protein_id\tproduct\n'+''.join(f'{v}\ttest\n' for v in ('A_g1.t1','A_g1.t2')))
    genes = root/'products.tsv'; genes.write_text('gene_id\tproduct\nA_g1\ttest\n')
    fs = status('functional.json',dict(status='ready',stage='functional_annotation',assembly_id='A',
        original_gff3_sha256=sha(identified/'A.gff3'),annotated_gff3_sha256=sha(coding),table_sha256=sha(table)))
    return dict(aid='A',identified=str(identified),genome=str(root/'genome.fa'),masked=str(masked),identity=str(identity),
        repeats=str(repeat),repeat_status=repeat_status,track_status=track_status,tracks=tracks,functional_enabled=True,
        functional=dict(gff=str(coding),table=str(table),genes=str(genes),status=fs))


class AnnotationReleaseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.request = release_fixture(self.root)

    def test_release_preserves_coding_and_maps_and_namespaces_repeats(self):
        out = self.root/'release'; result = package(self.request,out)
        self.assertEqual(result['genes'],1); self.assertEqual(result['coding_transcripts'],2)
        self.assertEqual(result['disabled_tracks'],['centromere'])
        self.assertEqual(result['sequence_counts'],dict(transcripts=2,trna=2))
        self.assertEqual(set(sequences(out/'A_transcripts.fa')),{'A_g1.t1','A_g1.t2'})
        self.assertEqual(sequences(out/'A_tRNA.fa'),{'A.trna':'A'*21,'A.trna.minus':'T'*17})
        for suffix in ('proteins','cds'):
            self.assertEqual(sha(out/f'A_{suffix}.fa'),sha(Path(self.request['identified'])/f'A_{suffix}.fa'))
        text = (out/'A.gff3').read_text()
        self.assertIn('ID=A.repeat.repeat1;Target=external_repeat',text)
        self.assertLess(text.index('ID=A.trna\n'), text.index('ID=A.trna.exon;'))
        self.assertIn('product=test',text)
        for name in ('gene_id_map.tsv','transcript_id_map.tsv','feature_id_map.tsv'):
            self.assertEqual(sha(out/'id_maps'/name),sha(Path(self.request['identified'])/name))
        for name,d in result['outputs'].items(): self.assertEqual(sha(out/name),d['sha256'])
        self.assertTrue((out/'tracks/A.infernal.gff3').is_file())
        self.assertFalse((out/'tracks/A.centromere.gff3').exists())
        checks = (out/'checksums.sha256').read_text()
        self.assertIn(sha(out/'release_status.json')+'  release_status.json', checks)
        self.assertEqual(len(gff_records(out/'A.gff3',{'chr1':1000})),sum(v['features'] for v in result['tracks'].values()))

    def test_no_functional_and_all_tracks_disabled(self):
        r = self.request; r['functional_enabled'] = False; r['functional'] = {}; r['tracks'] = []
        p = Path(r['track_status']); s = json.loads(p.read_text()); s['tracks'] = {t:dict(status='disabled') for t in s['tracks']}; p.write_text(json.dumps(s))
        state = package(r,self.root/'release')
        self.assertFalse(state['functional_enabled']); self.assertEqual(len(state['disabled_tracks']),4)
        self.assertFalse((self.root/'release/tables/functional_annotations.tsv').exists())
        self.assertEqual(state['sequence_counts'],dict(transcripts=2))
        self.assertFalse((self.root/'release/A_tRNA.fa').exists())

    def test_enabled_empty_trna_exports_empty_fasta(self):
        r=self.request;p=Path(r['tracks'][0]);p.write_text('##gff-version 3\n')
        s=Path(r['track_status']);state=json.loads(s.read_text())
        state['tracks']['trnascan'].update(features=0,gff3_sha256=sha(p));s.write_text(json.dumps(state))
        state=package(r,self.root/'empty')
        self.assertEqual(state['sequence_counts']['trna'],0)
        self.assertEqual((self.root/'empty/A_tRNA.fa').read_text(),'')
        self.assertIn('A_tRNA.fa',state['outputs'])

    def test_splicing_utr_both_strands_and_shared_exons(self):
        genome=self.root/'sequence.fa';genome.write_text('>chrX\naaCGTACCGGTT\n')
        targets={'plus':dict(seq='chrX',strand='+',start=0,end=12,exons=[(0,3),(8,12)]),
                 'minus':dict(seq='chrX',strand='-',start=0,end=12,exons=[(8,12),(0,3)])}
        destination=self.root/'tx.fa'
        extract_spliced(genome,[('transcripts',targets,destination)])
        # Exons are exported whole, including untranslated sequence; no CDS clipping.
        self.assertEqual(sequences(destination),dict(plus='AACGGTT',minus='AACCGTT'))
        rows=[(['chrX','native','tRNA','1','12','.','-','.',''],dict(ID='t')),
              (['chrX','native','exon','1','3','.','-','.',''],dict(ID='e1',Parent='t')),
              (['chrX','native','exon','9','12','.','-','.',''],dict(ID='e2',Parent='t'))]
        t=trna_models(rows)
        extract_spliced(genome,[('trna',t,destination)])
        self.assertEqual(sequences(destination),dict(t='AACCGTT'))
        # Native single-exon tRNAs may be represented by their feature span alone.
        t=trna_models(rows[:1]);extract_spliced(genome,[('trna',t,destination)])
        self.assertEqual(sequences(destination),dict(t='AACCGGTACGTT'))
        targets['plus']['exons']=[(0,4),(3,6)]
        with self.assertRaisesRegex(ValueError,'overlapping'):extract_spliced(genome,[('transcripts',targets,destination)])
        targets['plus']['exons']=[(0,3)];targets['plus']['strand']='.'
        with self.assertRaisesRegex(ValueError,'resolved strand'):extract_spliced(genome,[('transcripts',targets,destination)])

    def test_changed_track_and_missing_enabled_track_block(self):
        r = copy.deepcopy(self.request); r['tracks'] = r['tracks'][:-1]
        with self.assertRaisesRegex(ValueError,'Missing/extra selected'): package(r,self.root/'missing')
        Path(self.request['tracks'][0]).write_text('changed')
        with self.assertRaisesRegex(ValueError,'track changed'): package(self.request,self.root/'changed')

    def test_changed_coding_sequence_and_reference_block(self):
        p = self.root/'other.fa'; p.write_text((self.root/'genome.fa').read_text().replace('A','C'))
        r = copy.deepcopy(self.request); r['genome'] = str(p)
        with self.assertRaisesRegex(ValueError,'FASTA differs'): package(r,self.root/'reference_changed')
        p = Path(self.request['identified'])/'A_proteins.fa'; p.write_text(p.read_text()+'changed')
        with self.assertRaisesRegex(ValueError,'Identified output changed'): package(self.request,self.root/'changed')

    def test_collision_and_dangling_parent_block(self):
        p = Path(self.request['tracks'][0]); p.write_text(p.read_text().replace('A.trna','A_g1'))
        s = Path(self.request['track_status']); v = json.loads(s.read_text()); v['tracks']['trnascan']['gff3_sha256']=sha(p); s.write_text(json.dumps(v))
        with self.assertRaisesRegex(ValueError,'ID collision'): package(self.request,self.root/'collision')
        p.write_text('chr1\tx\texon\t1\t10\t.\t+\t.\tID=child;Parent=missing\n')
        with self.assertRaisesRegex(ValueError,'Dangling'): gff_records(p,{'chr1':1000})

    def test_parent_strand_bounds_and_cycles_block(self):
        p = self.root/'bad.gff3'
        for extra, pattern in [('chr1\tx\texon\t1\t10\t.\t-\t.\tID=e;Parent=p\n','child differs'),
                               ('chr1\tx\texon\t1\t50\t.\t+\t.\tID=e;Parent=p\n','child differs'),
                               ('chr1\tx\texon\t1\t10\t.\t+\t.\tID=e;Parent=e\n','Cyclic')]:
            p.write_text('chr1\tx\tgene\t1\t20\t.\t+\t.\tID=p\n'+extra)
            with self.assertRaisesRegex(ValueError,pattern): gff_records(p,{'chr1':1000})

    def test_overwrite_refused_and_checksums_stable(self):
        a=package(self.request,self.root/'one'); b=package(self.request,self.root/'two')
        self.assertEqual(a,b)
        with self.assertRaisesRegex(ValueError,'must be empty'): package(self.request,self.root/'one')


if __name__=='__main__': unittest.main()
