import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
import rna_alignment as alignment
from reference_repeats import reference as normalize


POLICY = dict(strandedness='auto', min_informative=100, dominant_fraction=.9, unstranded_tolerance=.1)


def qc(forward=100, reverse=0, status='ready'):
    return dict(assembly_id='A', readset_id='R', library_id='L', status=status,
                strand_forward=forward, strand_reverse=reverse)


class AlignmentContracts(unittest.TestCase):
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
            self.assertEqual(alignment.strand(qc(f, r), POLICY)['status'], 'blocked')

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
        result = alignment.report(plan, [], states, paths)
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['libraries'][0]['bams'], ['R.bam', 'R2.bam'])
        self.assertEqual(result['assemblies'][0]['bams'], ['R.bam', 'R2.bam'])

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
