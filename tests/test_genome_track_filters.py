import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from filter_genome_tracks import filter_track, trna_categories
from reference_repeats import sha


class TrackFilterTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.identity=self.root/'reference.json';self.identity.write_text(json.dumps(dict(reference_identity='test',contigs=[dict(id='chr1',length=500)])))

    def run_filter(self,track,lines,policy,aux=None):
        source=self.root/f'{track}.raw.gff3';source.write_text('##gff-version 3\n'+''.join(lines))
        status=self.root/'raw.json';status.write_text(json.dumps(dict(status='ready',assembly_id='A',track=track,reference_identity='test',gff3_sha256=sha(source),source_sha256='native',feature_types={})))
        result=filter_track('A',track,source,self.identity,status,policy,aux or [],self.root/'selected')
        self.assertEqual(result['gff3_sha256'],sha(self.root/f'selected/A.{track}.gff3'))
        return result

    def test_density_before_merge_preserves_dense_ends_but_rejects_internal(self):
        table=self.root/'windows.tsv'
        table.write_text('id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n'+''.join(f'chr1\t{end}\t{n}\t0\tTTAGGG\n' for end,n in [(100,10),(200,2),(300,10),(400,2),(500,10)]))
        lines=['chr1\ttidk\ttelomeric_repeat\t1\t500\t34\t+\t.\tID=raw\n']
        policy=dict(enabled=True,min_density=.5,max_end_distance=20,window=100,raw_min_count=2,motif='TTAGGG')
        s=self.run_filter('telomere',lines,policy,[table])
        self.assertEqual(s['features'],2);self.assertEqual(s['raw_features'],1)
        self.assertEqual(s['decision_counts'],dict(retained=2,rejected=3))
        self.assertEqual(s['reasons']['far_from_contig_end'],1)
        rows=[r.split('\t') for r in (self.root/'selected/A.telomere.gff3').read_text().splitlines()[1:]]
        self.assertEqual([(r[3],r[4]) for r in rows],[('1','100'),('401','500')])
        policy['min_density']=.8
        self.assertEqual(self.run_filter('telomere',lines,policy,[table])['features'],0)

    def test_partial_terminal_window_uses_actual_window_length(self):
        self.identity.write_text(json.dumps(dict(reference_identity='test',contigs=[dict(id='chr1',length=110)])))
        table=self.root/'windows.tsv';table.write_text('id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\nchr1\t100\t0\t0\tTTAGGG\nchr1\t110\t1\t0\tTTAGGG\n')
        s=self.run_filter('telomere',['chr1\ttidk\ttelomeric_repeat\t101\t110\t1\t+\t.\tID=raw\n'],dict(enabled=True,min_density=.5,max_end_distance=0,window=100,raw_min_count=1,motif='TTAGGG'),[table])
        self.assertEqual(s['features'],1)

    def test_native_trna_categories_and_children_follow_parent(self):
        table=self.root/'confidence.out';table.write_text('Sequence\t\t\t\t\t\nName\ttRNA#\tBegin\tEnd\tType\tNote\nchr1\t1\t1\t10\tAla\thigh confidence set\nchr1\t2\t30\t21\tAla\tsecondary filtered\n')
        lines=['chr1\ttRNAscan\ttRNA\t1\t10\t.\t+\t.\tID=t1\n','chr1\ttRNAscan\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1\n',
               'chr1\ttRNAscan\ttRNA\t21\t30\t.\t-\t.\tID=t2\n','chr1\ttRNAscan\texon\t21\t30\t.\t-\t.\tID=e2;Parent=t2\n']
        s=self.run_filter('trnascan',lines,dict(enabled=True),[table])
        self.assertEqual(s['features'],2);self.assertEqual(s['reasons']['parent_rejected'],1)
        self.assertEqual(self.run_filter('trnascan',lines,dict(enabled=False))['features'],4)
        table.write_text(table.read_text().replace('chr1\t1\t1\t10','chr1\t1\t2\t10'))
        with self.assertRaisesRegex(ValueError,'lacks tRNA coordinates'):self.run_filter('trnascan',lines,dict(enabled=True),[table])

    def test_trna_missing_native_note_is_an_error(self):
        table=self.root/'bad.out';table.write_text('chr1\t1\t1\t10\tAla\n')
        with self.assertRaisesRegex(ValueError,'Note column'):trna_categories(table)

    def test_rfam_uses_native_overlap_and_inclusion_without_new_evalue_cutoff(self):
        # All four native fmt-2 overlap categories, plus a below-inclusion hit.
        lines=[f'chr1\tInfernal\tncRNA\t{1+i*50}\t{40+i*50}\t30\t+\t.\tID=n{i};inc={inc};olp={olp};evalue=0.1\n' for i,(inc,olp) in enumerate([('!','^'),('!','='),('?','*'),('!','$'),('!','*')])]
        s=self.run_filter('infernal',lines,dict(enabled=True,remove_overlaps=True))
        self.assertEqual(s['features'],3)
        retained=(self.root/'selected/A.infernal.gff3').read_text()
        self.assertIn('ID=n0;',retained)
        self.assertIn('ID=n3;',retained)
        self.assertIn('ID=n4;',retained)
        self.assertEqual(s['reasons']['native_lower_ranked_overlap'],1)
        self.assertEqual(self.run_filter('infernal',lines,dict(enabled=True,remove_overlaps=False))['features'],4)
        for invalid in ('unknown','~'):
            with self.assertRaisesRegex(ValueError,'overlap annotation missing/unknown.*n0'):
                self.run_filter('infernal',[lines[0].replace('olp=^','olp='+invalid)],dict(enabled=True,remove_overlaps=True))
        with self.assertRaisesRegex(ValueError,'overlap annotation missing/unknown.*None'):
            self.run_filter('infernal',[lines[0].replace('olp=^;','')],dict(enabled=True,remove_overlaps=True))

    def test_centromere_identity_and_native_hor_confidence(self):
        lines=['chr1\tcentroAnno\tsatellite_DNA\t1\t50\t.\t.\t.\tID=m1;identity=0.9\n',
               'chr1\tcentroAnno\tsatellite_DNA\t51\t100\t.\t.\t.\tID=m2;identity=0.7\n',
               'chr1\tcentroAnno\trepeat_region\t101\t200\t.\t.\t.\tID=h1;identity=0.9;confidence=high\n',
               'chr1\tcentroAnno\trepeat_region\t201\t300\t.\t.\t.\tID=h2;identity=0.9;confidence=low\n']
        policy=dict(enabled=True,min_identity=.8,high_confidence_hor=True)
        s=self.run_filter('centromere',lines,policy);self.assertEqual(s['features'],2)
        self.assertEqual(s['reasons']['native_low_confidence_HOR'],1)
        with self.assertRaisesRegex(ValueError,'confidence missing'):self.run_filter('centromere',[lines[-1].replace('confidence=low','confidence=unknown')],policy)

    def test_empty_native_and_provenance_guard(self):
        self.assertEqual(self.run_filter('centromere',[],dict(enabled=True,min_identity=.8,high_confidence_hor=True))['features'],0)
        source=self.root/'centromere.raw.gff3';source.write_text('changed')
        with self.assertRaises(ValueError):filter_track('A','centromere',source,self.identity,self.root/'raw.json',dict(enabled=False),[],self.root/'selected')


if __name__=='__main__':unittest.main()
