import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
import rna_alignment as alignment
import report_rna_alignment as reporting
from reference_repeats import reference as normalize


POLICY = dict(strandedness='auto', min_informative=100, dominant_fraction=.9, unstranded_tolerance=.1)


def qc(forward=100, reverse=0, status='ready'):
    return dict(assembly_id='A', readset_id='R', library_id='L', status=status,
                strand_forward=forward, strand_reverse=reverse)


class AlignmentContracts(unittest.TestCase):
    def test_no_rna_report_link_resolves(self):
        import base64
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plan = root / 'plan.json'
            plan.write_text(json.dumps({'assemblies': [{'assembly_id': 'A', 'rna_readsets': []}], 'jobs': []}))
            encode = lambda value: base64.b64encode(json.dumps(value).encode()).decode()
            args = ['report_rna_alignment.py', '--plan', str(plan), '--qc', str(root / 'qc'),
                    '--strands', str(root / 'strands'), '--request', encode([]),
                    '--multiqc-request', encode(dict(paths=[], samples=[])), '--outdir', str(root)]
            with patch.object(sys, 'argv', args):
                reporting.main()
            self.assertIn('(multiqc/multiqc_report.html)', (root / 'README.md').read_text())
            self.assertIn('No RNA alignment was requested', (root / 'multiqc/multiqc_report.html').read_text())
            self.assertEqual(json.loads((root / 'rna_alignment_status.json').read_text())['status'], 'ready')

    def test_selected_read_bytes_and_mate_order(self):
        from prepare_rna import digest
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp) / 'R1.fq', Path(tmp) / 'R2.fq']
            for mate, path in enumerate(paths):
                path.write_text(f'fixture mate {mate}')
            spec = dict(layout='paired', sha256=[digest(path) for path in paths])
            alignment.selected_reads(spec, paths)
            with self.assertRaises(ValueError):
                alignment.selected_reads(spec, paths[::-1])
            paths[0].write_text('changed content')
            with self.assertRaises(ValueError):
                alignment.selected_reads(spec, paths)

    def test_strand_inference_and_ambiguity(self):
        for f, r, expected in ((100, 0, 'forward'), (90, 10, 'forward'),
                               (0, 100, 'reverse'), (10, 90, 'reverse'), (50, 50, 'unstranded')):
            self.assertEqual(alignment.strand(qc(f, r), POLICY)['strandedness'], expected)
        for f, r in ((1, 0), (75, 25), (0, 0)):
            self.assertEqual(alignment.strand(qc(f, r), POLICY)['status'], 'ready')
            self.assertEqual(alignment.strand(qc(f, r), {**POLICY, 'unresolved_policy': 'error'})['status'], 'blocked')

    def test_explicit_strand_does_not_hide_invalid_bam(self):
        explicit = {**POLICY, 'strandedness': 'reverse'}
        self.assertEqual(alignment.strand(qc(0, 0), explicit)['strandedness'], 'reverse')
        self.assertEqual(alignment.strand(qc(status='blocked'), explicit)['status'], 'blocked')

    def test_invalid_strand_thresholds(self):
        for override in ({'min_informative': 0}, {'dominant_fraction': .5}, {'unstranded_tolerance': .45}):
            with self.assertRaises(ValueError):
                alignment.strand(qc(), {**POLICY, **override})

    def test_library_conflict_is_isolated(self):
        plan = {'jobs': [dict(assembly_id='A', readset_id='R', library_id='L'),
                         dict(assembly_id='A', readset_id='R2', library_id='L')],
                'assemblies': [dict(assembly_id='A', rna_readsets=['R', 'R2']),
                               dict(assembly_id='B', rna_readsets=[])]}
        a = alignment.strand(qc(), POLICY)
        b = {**a, 'readset_id': 'R2', 'strandedness': 'reverse'}
        result = alignment.report(plan, [], [a, b], [{'assembly_id': 'A', 'readset_id': 'R'}])
        self.assertEqual([x['status'] for x in result['assemblies']], ['blocked', 'not_requested'])
        self.assertEqual(result['status'], 'incomplete')

    def test_no_rna_and_missing_bam(self):
        plan = {'jobs': [], 'assemblies': [dict(assembly_id='A', rna_readsets=[])]}
        self.assertEqual(alignment.report(plan, [], [], [])['status'], 'ready')
        plan['jobs'] = [dict(assembly_id='A', readset_id='R', library_id='L')]
        plan['assemblies'][0]['rna_readsets'] = ['R']
        self.assertEqual(alignment.report(plan, [qc()], [alignment.strand(qc(), POLICY)], [])['status'], 'incomplete')

    def test_bams_remain_separate_within_library(self):
        plan = {'jobs': [dict(assembly_id='A', readset_id=rid, library_id='L') for rid in ('R', 'R2')],
                'assemblies': [dict(assembly_id='A', rna_readsets=['R', 'R2'])]}
        states = [{**alignment.strand(qc(), POLICY), 'readset_id': rid} for rid in ('R', 'R2')]
        paths = [dict(assembly_id='A', readset_id=rid, bam=rid + '.bam', bam_index=rid + '.csi') for rid in ('R', 'R2')]
        result = alignment.report(plan, [{**qc(), 'readset_id': rid} for rid in ('R', 'R2')], states, paths)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['libraries'][0]['bams'], ['R.bam', 'R2.bam'])
        self.assertEqual(result['assemblies'][0]['bams'], ['R.bam', 'R2.bam'])

    def test_n2_fallback_preserves_evidence_and_strict_blocks(self):
        state = alignment.strand(qc(56916, 196953), POLICY)
        self.assertEqual(state['inference_outcome'], 'ambiguous')
        self.assertEqual(state['inferred_strandedness'], 'unresolved')
        self.assertEqual(state['informative'], 253869)
        self.assertAlmostEqual(state['reverse_fraction'], .7758, places=4)
        self.assertEqual(state['effective_strandedness'], 'unstranded')
        self.assertEqual(state['effective_setting_source'], 'conservative_fallback')
        self.assertIn('Intermediate', state['warning'])
        self.assertEqual(state['status'], 'ready')
        strict = alignment.strand(qc(56916, 196953), {**POLICY, 'unresolved_policy': 'error'})
        self.assertEqual(strict['status'], 'blocked')
        self.assertIn('strict strand policy', strict['reason'])

    def test_explicit_and_insufficient_provenance(self):
        for setting in ('forward', 'reverse', 'unstranded'):
            state = alignment.strand(qc(56916, 196953), {**POLICY, 'strandedness': setting})
            self.assertEqual(state['effective_setting_source'], 'explicit')
            self.assertEqual(state['effective_strandedness'], setting)
            self.assertEqual(state['inferred_strandedness'], 'unresolved')
            self.assertEqual(state['strand_reverse'], 196953)
        state = alignment.strand(qc(2, 1), POLICY)
        self.assertEqual(state['inference_outcome'], 'insufficient')
        self.assertIn('Only 3', state['warning'])
        for policy in ('unstranded', 'error'):
            state = alignment.strand(qc(status='blocked'), {**POLICY, 'unresolved_policy': policy})
            self.assertEqual(state['status'], 'blocked')
            self.assertEqual(state['inference_outcome'], 'not_evaluated')
        for override in ({'strandedness': 'wrong'}, {'unresolved_policy': 'wrong'}):
            with self.assertRaises(ValueError): alignment.strand(qc(), {**POLICY, **override})

    def library_report(self, policies, counts):
        ids = ['R', 'R2']
        plan = dict(jobs=[dict(assembly_id='A', readset_id=rid, library_id='L') for rid in ids],
                    assemblies=[dict(assembly_id='A', rna_readsets=ids)])
        statistics = [{**qc(*count), 'readset_id': rid} for rid, count in zip(ids, counts)]
        states = [alignment.strand(row, policy) for row, policy in zip(statistics, policies)]
        paths = [dict(assembly_id='A', readset_id=rid, bam=rid+'.bam', bam_index=rid+'.csi') for rid in ids]
        return alignment.report(plan, statistics, states, paths)

    def test_lane_disagreement_and_unresolved_library_fallback(self):
        for counts in ([(100, 0), (0, 100)], [(100, 0), (56916, 196953)]):
            result = self.library_report([POLICY, POLICY], counts)
            self.assertEqual(result['status'], 'ready')
            self.assertEqual(result['libraries'][0]['effective_strandedness'], 'unstranded')
            self.assertTrue(result['has_warnings'])
            self.assertEqual(len(result['warnings']), 2)
            self.assertTrue(all(r['strand']['effective_strandedness'] == 'unstranded' for r in result['readsets']))
            self.assertEqual(result['readsets'][0]['strand']['inferred_strandedness'], 'forward')
            self.assertIn('conservative_fallback', alignment.markdown(result))
        strict = {**POLICY, 'unresolved_policy': 'error'}
        self.assertEqual(self.library_report([strict, strict], [(100, 0), (0, 100)])['status'], 'incomplete')
        explicit = [{**POLICY, 'strandedness': setting} for setting in ('forward', 'reverse')]
        self.assertEqual(self.library_report(explicit, [(100, 0), (0, 100)])['status'], 'incomplete')

    def test_missing_or_invalid_alignment_cannot_be_rescued_by_fallback(self):
        state = alignment.strand(qc(56916, 196953), POLICY)
        plan = dict(jobs=[dict(assembly_id='A', readset_id='R', library_id='L')],
                    assemblies=[dict(assembly_id='A', rna_readsets=['R'])])
        paths = [dict(assembly_id='A', readset_id='R', bam='r.bam', bam_index='r.csi')]
        for statistics, strands in (([], [state]), ([qc(status='blocked')], [state]), ([qc()], [])):
            result = alignment.report(plan, statistics, strands, paths)
            self.assertEqual(result['status'], 'incomplete')
            self.assertIn('blocked', alignment.markdown(result))

    def test_reference_content_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source.fa'
            source.write_text('>chr1\nACGTACGT\n')
            normalize(source, 'A', 0, root / 'A')
            identity = alignment.load(root / 'A/A.contigs.json')
            fasta = root / 'A/A.filtered.fasta'
            self.assertEqual(alignment.reference(fasta, identity)['reference_identity'], identity['reference_identity'])
            fasta.write_text('>chr1\nACGTACGA\n')
            with self.assertRaises(ValueError):
                alignment.reference(fasta, identity)

    def test_plan_keeps_reference_specific_jobs_and_rejects_changed_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'source.fa'
            reads = root / 'reads.fq'
            fasta.write_text('>chr1\nACGTACGT\n')
            reads.write_text('@r\nACGT\n+\nIIII\n')
            normalize(fasta, 'A', 0, root / 'A')
            normalize(fasta, 'B', 0, root / 'B')
            row = dict(readset_id='R', library_id='L', sample_id=None, layout='single',
                       status='ready', read_paths=[str(reads)], read_units=1, sha256=['fixture'])
            manifest = dict(status='ready', rna_readsets=[{**row, 'strandedness': 'auto'}],
                            assemblies=[dict(assembly_id=x, rna_readsets=['R']) for x in ('A', 'B')])
            prepared = dict(status='ready', readsets=[row])
            locations = dict(reference_dir=str(root), strand_policy=POLICY)
            planned = alignment.plan(manifest, prepared, locations)
            self.assertEqual(len(planned['jobs']), 2)
            self.assertEqual(len({x['reference_id'] for x in planned['references']}), 1)
            self.assertNotIn('strandedness', planned['jobs'][0]['spec'])
            row['library_id'] = 'wrong'
            with self.assertRaises(ValueError):
                alignment.plan(manifest, prepared, locations)


class BamStatsContracts(unittest.TestCase):
    def run_stats(self, reads, expected=1, references=('chr1',), lengths=(1000,)):
        class Bam:
            def __enter__(self):
                return self
            def __exit__(self, *args):
                pass
            def fetch(self, until_eof):
                return iter(reads)
        bam = Bam()
        bam.references, bam.lengths = references, lengths
        bam.header = SimpleNamespace(to_dict=lambda: {'RG': [dict(ID='R', SM='L', LB='L')]})
        identity = dict(reference_identity='id', contigs=[dict(id='chr1', length=1000)])
        spec = dict(assembly_id='A', readset_id='R', library_id='L', layout='single',
                    read_units=expected, reference_identity='id', min_mapq=20)
        with tempfile.TemporaryDirectory() as tmp, patch.dict(sys.modules, {'pysam': SimpleNamespace(AlignmentFile=lambda *a: bam)}):
            alignment.stats('fixture.bam', identity, spec, tmp)
            return alignment.load(Path(tmp) / 'alignment.json')

    @staticmethod
    def read(reverse=False, xs='+', position=10, nh=1):
        tags = {'XS': xs, 'NH': nh, 'RG': 'R'}
        return SimpleNamespace(reference_id=0, reference_start=position, reference_end=position + 150,
                               reference_name='chr1', is_secondary=False, is_supplementary=False,
                               is_paired=False, is_read1=False, is_unmapped=False, is_reverse=reverse,
                               is_proper_pair=False, mapping_quality=60, cigartuples=[(0, 25), (3, 100), (0, 25)],
                               has_tag=lambda tag: tag in tags, get_tag=lambda tag: tags[tag])

    def test_orientation_and_unique_splice_support(self):
        for reverse, xs, expected in ((False, '+', 'forward'), (True, '-', 'forward'),
                                      (True, '+', 'reverse'), (False, '-', 'reverse')):
            result = self.run_stats([self.read(reverse, xs)])
            self.assertEqual(result['status'], 'ready')
            self.assertEqual(result['strand_' + expected], 1)
            self.assertEqual(result['junctions'], 1)
        result = self.run_stats([self.read(nh=2)])
        self.assertEqual(result['strand_forward'], 0)
        self.assertEqual(result['junctions'], 0)

    def test_xs_availability_allows_ambiguous_records(self):
        result = self.run_stats([self.read(xs='.'), self.read(xs='+', position=20)], expected=2)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['spliced_xs_reads'], 1)
        result = self.run_stats([self.read(xs='.')])
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['spliced_xs_reads'], 0)

    def test_dictionary_sort_and_counts(self):
        self.assertEqual(self.run_stats([self.read()], references=('wrong',))['status'], 'blocked')
        self.assertEqual(self.run_stats([self.read(position=20), self.read(position=10)], expected=2)['status'], 'blocked')
        self.assertEqual(self.run_stats([self.read()], expected=2)['status'], 'blocked')
        self.assertEqual(self.run_stats([])['status'], 'blocked')

    def test_wrong_readgroup_and_mate_flags(self):
        read = self.read()
        read.get_tag = lambda tag: 'wrong' if tag == 'RG' else (1 if tag == 'NH' else '+')
        self.assertEqual(self.run_stats([read])['status'], 'blocked')
        read = self.read()
        read.is_paired = True
        self.assertEqual(self.run_stats([read])['status'], 'blocked')


if __name__ == '__main__':
    unittest.main()
