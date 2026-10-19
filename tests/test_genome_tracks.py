import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from genome_tracks import TRACKS, centromeres, normalize, telomeres
from report_genome_tracks import report
from validate_genome_tracks_cluster import task_process, fixtures, check_telomere_ends
from reference_repeats import fasta


class GenomeTracksTests(unittest.TestCase):
    def test_satellite_motif_can_join_both_telomere_arrays(self):
        fixtures(self.root)
        genome=self.root/'fixtures/B.fa';tsv=self.root/'windows.tsv'
        motif='TTAGG';reverse='CCTAA';window=1000
        with tsv.open('w') as handle:
            handle.write('id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n')
            for name,seq in fasta(genome):
                for start in range(0,len(seq),window):
                    part=seq[start:start+window]
                    handle.write(f'{name}\t{start+len(part)}\t{part.count(motif)}\t{part.count(reverse)}\t{motif}\n')
        output=self.root/'B.telomere.gff3'
        telomeres(tsv,genome,'B',motif,window,2,1000,output)
        rows=[s for s in output.read_text().splitlines() if not s.startswith('#')]
        self.assertEqual(len(rows),1)
        self.assertIn('position=spanning',rows[0])
        check_telomere_ends(self.root,'B',self.root/'fixtures/B.json',motif)
        output.write_text(output.read_text().replace('\t1\t','\t2\t',1))
        with self.assertRaisesRegex(AssertionError,'missing left'):check_telomere_ends(self.root,'B',self.root/'fixtures/B.json',motif)

    def test_exact_process_match_excludes_infernal_merge(self):
        tasks=[{'name':f'GENOME_TRACKS:INFERNAL_PARALLEL:INFERNAL ({aid}:chunk_{chunk:03d})'}
               for aid in 'AB' for chunk in range(2)]
        tasks += [{'name':f'GENOME_TRACKS:INFERNAL_PARALLEL:MERGE_INFERNAL ({aid})'} for aid in 'AB']
        self.assertEqual(sum(task_process(t)=='INFERNAL' for t in tasks),4)
        self.assertEqual(sum(task_process(t)=='MERGE_INFERNAL' for t in tasks),2)
        self.assertEqual(task_process({'name':'GENOME_TRACKS:CENTROMERE_TO_GFF3 (A)'}),'CENTROMERE_TO_GFF3')

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.identity = self.root/'identity.json'
        self.identity.write_text(json.dumps(dict(reference_identity='test',contigs=[dict(id='chr1',length=25)])))

    def test_tidk_end_coordinates_partial_window_and_terminal(self):
        genome = self.root/'genome.fa'; genome.write_text('>chr1\n'+'A'*25+'\n')
        tsv = self.root/'tidk.tsv'
        tsv.write_text('id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\nchr1\t10\t2\t0\tTTAGGG\nchr1\t20\t0\t0\tTTAGGG\nchr1\t25\t0\t3\tTTAGGG\n')
        out = self.root/'telomere.gff'
        telomeres(tsv,genome,'A','TTAGGG',10,2,2,out)
        rows=[s.split('\t') for s in out.read_text().splitlines() if not s.startswith('#')]
        self.assertEqual([(r[3],r[4]) for r in rows],[('1','10'),('21','25')])
        self.assertIn('position=5_prime',rows[0][8]);self.assertIn('position=3_prime',rows[1][8])
        tsv.write_text(tsv.read_text().replace('chr1\t10','chr1\t0'))
        with self.assertRaisesRegex(ValueError,'window end'):telomeres(tsv,genome,'A','TTAGGG',10,2,2,out)

    def test_centromere_native_inclusive_end_headers_and_bad_length(self):
        mono=self.root/'mono.csv';hor=self.root/'hor.csv';out=self.root/'centro.gff'
        mono.write_text('sequence name,mononer name,start position,end position,estimated identity,length\nchr1,Pos:0:25_1,0,9,0.95,10\n')
        hor.write_text('chr1,Pos:0:25:1_2,10,24,0.90,15,high\n')
        centromeres(mono,hor,'A',out)
        rows=[s.split('\t') for s in out.read_text().splitlines() if not s.startswith('#')]
        self.assertEqual([(r[3],r[4]) for r in rows],[('1','10'),('11','25')])
        self.assertIn('confidence=high',rows[1][8])
        mono.write_text('chr1,x,0,9,0.95,9\n')
        with self.assertRaisesRegex(ValueError,'coordinates/length'):centromeres(mono,hor,'A',out)

    def test_gff_parents_bounds_duplicates_and_rfam_family(self):
        native=self.root/'native.gff';out=self.root/'normalized.gff'
        native.write_text('chr1\tInfernal\tTest_RNA\t1\t25\t10\t-\t.\tevalue=0.001\n')
        # Status sidecar is emitted beside the process output in its working directory.
        with patch('genome_tracks.dump') as dump:
            state=normalize('A','infernal',native,self.identity,out)
        self.assertEqual(state['features'],1);self.assertIn('\tncRNA\t',out.read_text());self.assertIn('Name=Test_RNA',out.read_text())
        native.write_text('chr1\ttest\ttRNA\t1\t26\t.\t+\t.\tID=x\n')
        with self.assertRaisesRegex(ValueError,'outside reference'):normalize('A','trnascan',native,self.identity,out)
        native.write_text('chr1\ttest\texon\t1\t25\t.\t+\t.\tID=x;Parent=missing\n')
        with self.assertRaisesRegex(ValueError,'parents'):normalize('A','trnascan',native,self.identity,out)
        native.write_text(('chr1\ttest\ttRNA\t1\t25\t.\t+\t.\tID=x\n')*2)
        with self.assertRaisesRegex(ValueError,'Duplicate'):normalize('A','trnascan',native,self.identity,out)

    def test_disabled_empty_and_missing_results(self):
        settings={'run_'+t:False for t in TRACKS}
        state=report('A',settings,self.identity,[],self.root/'report')
        self.assertEqual(state['tracks']['infernal']['status'],'disabled')
        native=self.root/'empty.gff';native.write_text('##gff-version 3\n')
        with patch('genome_tracks.dump'):state=normalize('A','infernal',native,self.identity,self.root/'empty.out')
        self.assertEqual(state['features'],0)
        settings['run_infernal']=True
        with self.assertRaisesRegex(ValueError,'Missing/duplicate'):report('A',settings,self.identity,[],self.root/'report')
        status=self.root/'status.json';status.write_text(json.dumps(state))
        report('A',settings,self.identity,[status],self.root/'report')
        state['reference_identity']='wrong';status.write_text(json.dumps(state))
        with self.assertRaisesRegex(ValueError,'identity mismatch'):report('A',settings,self.identity,[status],self.root/'report')


if __name__ == '__main__': unittest.main()
