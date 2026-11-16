from collections import Counter
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from reference_repeats import sha
from report_utr_review import junction_counts, junction_flags, mirror_matches, read_table, review
from unified_report import utr_review_section


class UtrReviewReportTests(unittest.TestCase):
    def test_exact_junction_counts_respect_strand_anchor_and_unique_mapping(self):
        class Read(SimpleNamespace):
            def has_tag(self, tag):
                return tag == 'NH'

            def get_tag(self, tag):
                return self.nh

        def read(reverse=False, nh=1, anchor=10):
            return Read(is_reverse=reverse, nh=nh, cigartuples=[(0,10),(3,40),(0,anchor)],
                        reference_start=90, mapping_quality=60, is_unmapped=False,
                        is_secondary=False, is_supplementary=False, is_qcfail=False, is_paired=False)

        bam = SimpleNamespace(fetch=lambda *args: [read(), read(reverse=True), read(nh=2), read(anchor=7)])
        policy = dict(anchor_bp=8, min_mapq=20)
        self.assertEqual(junction_counts(bam, 'chr', (100,140), '-', 'reverse', policy),
                         Counter(proposed=1, opposite=1))
        self.assertEqual(junction_counts(bam, 'chr', (100,140), '-', 'unstranded', policy),
                         Counter(unstranded=2))

    def test_same_strand_junction_counts_distinguish_zero_and_distributed_support(self):
        stranded = {'L1', 'L2'}
        self.assertEqual(junction_flags({'L1': Counter(opposite=8), 'L2': Counter(unstranded=6)}, stranded, 3),
                         {'new_junction_supported_only_by_unstranded_reads'})
        self.assertEqual(junction_flags({'L1': Counter(proposed=2), 'L2': Counter(proposed=2)}, stranded, 3),
                         {'weak_stranded_junction_support_across_libraries'})

    def test_mirror_requires_reciprocal_exonic_overlap_and_same_splice_chain(self):
        rows = {'t': dict(contig='chr', strand='+', evidence_id='sense')}
        models = {('chr', '+', 'sense'): [(10,20), (90,180)],
                  ('chr', '-', 'exact'): [(10,20), (90,180)],
                  ('chr', '-', 'near'): [(9,20), (90,179)],
                  ('chr', '-', 'different_junction'): [(10,21), (90,180)],
                  ('chr', '-', 'shorter'): [(10,20), (90,150)]}
        matches, missing = mirror_matches(rows, models)
        self.assertEqual(matches['t'], [('exact','exact'), ('near','near')])
        self.assertFalse(missing)

    def test_review_links_release_ids_and_preserves_strand_across_libraries(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            release = root/'release'
            (release/'id_maps').mkdir(parents=True)
            (release/'provenance').mkdir()

            def put(path, value):
                path.write_text(value, encoding='utf-8')
                return path

            def write_json(path, value):
                return put(path, json.dumps(value))

            def gff(feature, start, end, strand, attrs):
                return f'chr\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n'

            def model(tid, gene, start, end, strand):
                return (gff('mRNA',start,end,strand,f'ID={tid};Parent={gene}') +
                        gff('exon',start,end,strand,f'Parent={tid}') +
                        gff('CDS',start,end,strand,f'Parent={tid}'))

            original = put(root/'original.gff3', model('t','g',101,120,'-') +
                           model('nearby','neighbor',56,60,'+'))
            screened = put(root/'screened.gff3', original.read_text() +
                           gff('exon',51,70,'-','Parent=t'))
            clustered = put(root/'clustered.gtf', gff('exon',51,70,'-','transcript_id "evidence"'))
            txmap = put(release/'id_maps/transcript_id_map.tsv',
                        'transcript_id\tdecision\tfinal_transcript_id\nt\tretained\tfinal.t1\n')
            gmap = put(release/'id_maps/gene_id_map.tsv', 'gene_id\tfinal_gene_id\nneighbor\tfinal.neighbor\n')
            ref = write_json(release/'provenance/reference.json', dict(contigs=[dict(id='chr',length=200)]))
            rel_status = write_json(release/'release_status.json', dict(status='ready',reference_identity='ref',
                outputs={p.relative_to(release).as_posix(): dict(sha256=sha(p)) for p in (txmap,gmap,ref)}))
            libs = [('U','unstranded'), ('L1','reverse'), ('L2','forward')]
            status = write_json(root/'utr.json', dict(reference_identity='ref',original_model_sha256=sha(original),
                gff3_sha256=sha(screened),utr_support_policy=dict(min_breadth=.8,min_junction_alignments=3,anchor_bp=8,min_mapq=20),
                readsets=[dict(library_id=k,readset_id=k,effective_strandedness=v) for k,v in libs]))
            per_library = [dict(library_id=k,covered_fraction=.9,all_new_junctions_supported=k=='U',
                all_boundary_links_supported=True,min_new_junction_alignments=n)
                for (k,_),n in zip(libs,(3,2,1))]
            support_fields = ['transcript_id','gene_id','contig','strand','selected_model','new_junctions',
                'opposite_strand_cds_overlap_bp','added_5prime_spliced_bp','added_3prime_spliced_bp',
                'added_blocks_0based_halfopen','evidence_id','library_support_json']
            support = put(root/'support.tsv', '\t'.join(support_fields)+'\n'+
                '\t'.join(['t','g','chr','-','new','1','5','0','20','[[50,70]]','evidence',json.dumps(per_library)])+'\n')

            class Bam:
                references, lengths = ['chr'], [200]

                def __init__(self, filename, *args, **kwargs):
                    self.filename = filename

                def __enter__(self):
                    return self

                def __exit__(self, *args):
                    return False

            def junction(bam, *args):
                return {'U': Counter(unstranded=3), 'L1': Counter(proposed=2,opposite=1),
                        'L2': Counter(proposed=1)}[bam.filename]

            out = root/'UTR'
            with patch.dict(sys.modules, pysam=SimpleNamespace(AlignmentFile=Bam)), \
                    patch('report_utr_review.junction_counts', side_effect=junction), \
                    patch('report_utr_review.overlap_depth', return_value=(5,10)):
                review(original, screened, support, clustered, release, status,
                       [k for k,_ in libs], [k+'.csi' for k,_ in libs], out)
            flags = read_table(out/'utr_review.tsv')
            self.assertEqual(flags[0]['release_transcript_id'], 'final.t1')
            self.assertIn('weak_stranded_junction_support_across_libraries', flags[0]['flags'])
            junctions = {r['library']: r for r in read_table(out/'utr_junction_support.tsv')}
            self.assertEqual(junctions['L1']['opposite_strand_alignments'], '1')
            self.assertEqual(junctions['U']['proposed_strand_alignments'], '.')
            overlap = read_table(out/'utr_cds_overlap_coverage.tsv')
            self.assertEqual((overlap[0]['opposite_release_gene_id'], overlap[0]['start_0based']), ('final.neighbor','55'))
            data = root/'report_data'
            data.mkdir()
            section = utr_review_section(out, data, 'assembly', 'ref', sha(rel_status), {})
            self.assertIn('assembly.utr_review.tsv','\n'.join(section))
            self.assertEqual(read_table(data/'assembly.utr_review.tsv'), flags)
