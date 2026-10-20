import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from test_final_structural import accepted
from package_annotation import package, gff_records
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
        text = {'trnascan':'chr1\ttRNAscan\texon\t10\t30\t.\t+\t.\tID=A.trna.exon;Parent=A.trna\nchr1\ttRNAscan\ttRNA\t10\t30\t.\t+\t.\tID=A.trna\n',
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
