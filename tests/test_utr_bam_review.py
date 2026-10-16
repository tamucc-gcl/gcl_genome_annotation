from collections import Counter
from types import SimpleNamespace
import unittest

from review_utr_bams import cigar_parts, depth_metrics, library_strand, inspect_bam


class UtrBamReview(unittest.TestCase):
    def test_library_orientation_for_both_mates_and_protocols(self):
        for protocol in ('forward','reverse'):
            for target in ('+','-'):
                read1_reverse=(target=='-') != (protocol=='reverse')
                for mate, reverse in [(1,read1_reverse),(2,not read1_reverse)]:
                    read=SimpleNamespace(is_paired=True,is_reverse=reverse,is_read1=mate==1,is_read2=mate==2)
                    self.assertEqual(library_strand(read,protocol),target)
                    self.assertIsNone(library_strand(read,'unstranded'))

    def test_cigar_deletions_are_not_coverage_or_introns(self):
        read=SimpleNamespace(reference_start=100,cigartuples=[(4,5),(0,10),(3,90),(0,8),(2,2),(0,10)])
        matches,gaps=cigar_parts(read)
        self.assertEqual(matches,[(100,110),(200,208),(210,220)])
        self.assertEqual(gaps,[(110,200,8)])

    def test_zero_coverage_and_intron_gaps_are_handled_separately(self):
        events=Counter({10:2,15:-1,20:-1,30:1,35:-1})
        m=depth_metrics(events,[(10,25),(30,40)])
        self.assertEqual(m['bases'],25)
        self.assertEqual(m['mean_depth'],20/25)
        self.assertEqual(m['covered_fraction'],15/25)
        self.assertEqual(m['longest_uncovered_run_bp'],5)
        self.assertEqual(depth_metrics({},[(0,20)])['longest_uncovered_run_bp'],20)
        self.assertIsNone(depth_metrics({},[])['mean_depth'])

    def test_distinct_fragments_anchor_filter_and_strand_counts(self):
        def read(name, reverse=True, anchor=10, nh=1):
            tags={'NH':nh,'XS':'+'}
            r=SimpleNamespace(query_name=name,reference_start=10,cigartuples=[(0,10),(3,10),(0,anchor)],
                is_paired=False,is_reverse=reverse,is_unmapped=False,is_secondary=False,
                is_supplementary=False,is_qcfail=False,mapping_quality=60)
            r.has_tag=lambda key:key in tags
            r.get_tag=lambda key:tags[key]
            return r
        class Bam:
            def fetch(self,*args):return iter([read('pair'),read('pair'),read('other',False),read('short',anchor=5),read('multi',nh=2)])
        case=dict(accepted=dict(seq='chr1',strand='+',exons=[(10,20),(30,40)]),
                  original=dict(exons=[(10,20)]),accepted_extension=[(30,40)])
        result=inspect_bam(Bam(),case,dict(readset_id='R',strandedness='reverse'),20,8)
        gap=result['expected_junctions'][0]
        self.assertTrue(gap['new_in_refined_model'])
        self.assertEqual(gap['supporting_fragments']['all_strands'],2)
        self.assertEqual(gap['supporting_fragments']['matching_library_strand'],1)
        self.assertEqual(gap['supporting_fragments']['opposite_library_strand'],1)
        self.assertEqual(result['extension_depth']['matching_library_strand']['mean_depth'],2.5)
