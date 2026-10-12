import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
import braker_prediction as braker
from report_braker import summarize
from reference_repeats import reference


def inputs():
    preflight = dict(status='ready', assemblies=[dict(assembly_id='A', prediction_mode='rna_protein',
                    rna_readsets=['R2', 'R1'], settings={'genetic_code': 1})])
    repeats = dict(status='ready', assemblies=[dict(assembly_id='A', status='ready',
                   repeats=dict(reference_identity='ref', masked_sha256='genome'))])
    proteins = dict(status='ready', assemblies=[dict(assembly_id='A', status='ready',
                    protein_count=2, proteins_sha256='prot', source_ids=['db', 'spike'])])
    rows = [dict(assembly_id='A', readset_id=rid, library_id=rid, sample_id=rid, layout='paired',
                 status='ready', reference_identity='ref', spliced_primary_reads=100,
                 strand=dict(status='ready', strandedness='reverse')) for rid in ('R1', 'R2')]
    preflight['rna_readsets'] = [dict(readset_id=r['readset_id'], library_id=r['library_id'],
                                    sample_id=r['sample_id'], layout=r['layout'], strandedness='auto') for r in rows]
    libs = [dict(assembly_id='A', library_id=r['library_id'], status='ready', readsets=[r['readset_id']],
                 bams=[r['readset_id']+'.bam'], bam_indices=[r['readset_id']+'.bam.csi']) for r in rows]
    return preflight, repeats, proteins, dict(status='ready', readsets=rows, libraries=libs)


class BrakerHandoff(unittest.TestCase):
    def test_separate_bams_keep_stable_order_and_reverse_metadata(self):
        job = braker.plan(*inputs())['jobs'][0]
        self.assertEqual([r['readset_id'] for r in job['readsets']], ['R1', 'R2'])
        self.assertEqual([r['bam'] for r in job['readsets']], ['R1.bam', 'R2.bam'])
        self.assertEqual([r['strandedness'] for r in job['readsets']], ['reverse', 'reverse'])
        args = braker.command(job, 'genome.fa', 'proteins.fa', ['R1.bam', 'R2.bam'], 4, 8, 'braker', 'augustus')
        bam_arg = next(a for a in args if a.startswith('--bam='))
        self.assertEqual(bam_arg.count(','), 1)
        self.assertIn('--rounds=8', args)
        self.assertFalse(any(a.startswith(('--stranded', '--UTR', '--busco_lineage')) for a in args))
        self.assertEqual(args[-1], '--verbosity=3')

    def test_incomplete_reference_and_wrong_assembly_identity_block(self):
        for index in range(4):
            data = inputs()
            data[index]['status'] = 'incomplete'
            with self.assertRaises(ValueError): braker.plan(*data)
        data = inputs(); data[3]['readsets'][0]['reference_identity'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'identities differ'): braker.plan(*data)

    def test_unresolved_strand_empty_splice_and_unsupported_mode_block(self):
        for change in ('strand', 'splice', 'mode'):
            data = inputs()
            if change == 'strand': data[3]['readsets'][0]['strand']['strandedness'] = 'auto'
            elif change == 'splice': data[3]['readsets'][0]['spliced_primary_reads'] = 0
            else: data[0]['assemblies'][0]['prediction_mode'] = 'protein'
            with self.assertRaises(ValueError): braker.plan(*data)

    def test_library_bam_count_and_conflicting_lane_strands_block(self):
        data = inputs(); data[3]['libraries'][0]['bam_indices'] = []
        with self.assertRaises(ValueError): braker.plan(*data)
        data = inputs()
        for r in data[3]['readsets']: r['library_id'] = 'L'
        for r in data[0]['rna_readsets']: r['library_id'] = 'L'
        data[3]['readsets'][1]['strand']['strandedness'] = 'forward'
        with self.assertRaises(ValueError): braker.plan(*data)

    def test_software_check_does_not_request_prediction(self):
        args = braker.command(braker.plan(*inputs())['jobs'][0], 'g', 'p', ['b'], 1, 8, 'w', 'a', check=True)
        self.assertIn('--checkSoftware', args)
        with self.assertRaises(ValueError): braker.command({'mode': 'rna_protein'}, 'g', 'p', [], 1, 8, 'w', 'a')
        with self.assertRaises(ValueError): braker.command({'mode': 'rna_protein'}, 'g', 'p', ['b'], 64, 8, 'w', 'a')

    def test_nonstandard_genetic_code_and_reused_bam_block(self):
        data = inputs(); data[0]['assemblies'][0]['settings']['genetic_code'] = 6
        with self.assertRaisesRegex(ValueError, 'standard genetic code'): braker.plan(*data)
        data = inputs(); data[3]['libraries'][1]['bams'] = ['R1.bam']
        with self.assertRaisesRegex(ValueError, 'same BAM'): braker.plan(*data)

    def test_current_sheet_must_match_accepted_bam_metadata_and_strand(self):
        data = inputs(); data[0]['rna_readsets'][0]['sample_id'] = 'other'
        with self.assertRaisesRegex(ValueError, 'RNA sheet differs'): braker.plan(*data)
        data = inputs(); data[0]['rna_readsets'][0]['strandedness'] = 'forward'
        with self.assertRaisesRegex(ValueError, 'strand setting differs'): braker.plan(*data)


class PredictionIntegrity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.identity = dict(reference_identity='ref', contigs=[dict(id='chr1', length=100)])
        self.spec = dict(assembly_id='A', reference_identity='ref')
        self.gff = self.root/'genes.gff3'; self.aa = self.root/'aa.fa'; self.cds = self.root/'cds.fa'
        self.gff.write_text('##gff-version 3\nchr1\tBRAKER3\tgene\t1\t9\t.\t+\t.\tID=g1\n'
            'chr1\tBRAKER3\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n'
            'chr1\tBRAKER3\tCDS\t1\t9\t.\t+\t0\tID=g1.t1.CDS;Parent=g1.t1\n')
        self.aa.write_text('>g1.t1\nMFK\n'); self.cds.write_text('>g1.t1\nATGTTTAAA\n')

    def validate(self):
        return braker.validate_outputs(self.spec, self.identity, self.gff, self.aa, self.cds)

    def test_nonempty_gff_and_matching_sequence_outputs(self):
        self.assertEqual(self.validate()['coding_transcripts'], 1)

    def test_out_of_bounds_and_missing_parent_block(self):
        original = self.gff.read_text()
        for text in (original.replace('\t9\t', '\t101\t'), original.replace('Parent=g1\n', 'Parent=missing\n')):
            self.gff.write_text(text)
            with self.assertRaises(ValueError): self.validate()

    def test_missing_or_duplicate_predicted_sequences_block(self):
        for text in ('>wrong\nMFK\n', '>g1.t1\nMFK\n>g1.t1\nMFK\n', ''):
            self.aa.write_text(text)
            with self.assertRaises(ValueError): self.validate()

    def test_report_distinguishes_runtime_check_from_predictions(self):
        plan = braker.plan(*inputs())
        check = summarize(plan, [], [dict(assembly_id='A')], True, self.root)
        self.assertEqual(check['status'], 'software_checked')
        self.assertNotIn('(braker/A/A.gff3)', (self.root/'README.md').read_text())
        missing = summarize(plan, [], [dict(assembly_id='A')], False, self.root)
        self.assertEqual(missing['status'], 'incomplete')
        ready = summarize(plan, [self.validate()], [dict(assembly_id='A')], False, self.root)
        self.assertEqual(ready['status'], 'ready')
        self.assertIn('(braker/A/A.gff3)', (self.root/'README.md').read_text())


class BamHandoffIntegrity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        raw = self.root/'raw.fa'; raw.write_text('>chr1\nACGTACGT\n')
        reference(raw, 'A', 0, self.root/'ref')
        self.genome = self.root/'ref/A.filtered.fasta'
        self.identity = braker.load(self.root/'ref/A.contigs.json')
        self.protein = self.root/'proteins.fa'; self.protein.write_text('>p\nMFK\n')
        self.index = self.root/'reads.bam.csi'; self.index.write_text('mock index')
        self.spec = dict(assembly_id='A', reference_identity=self.identity['reference_identity'],
                         masked_sha256=braker.sha(self.genome), proteins_sha256=braker.sha(self.protein),
                         readsets=[dict(readset_id='R', library_id='L', sample_id='S')])
        self.header = {'HD': {'SO': 'coordinate'}, 'RG': [{'ID': 'R', 'LB': 'L', 'SM': 'S'}]}
        self.bam = MagicMock()
        self.bam.references = ('chr1',); self.bam.lengths = (8,)
        self.bam.header.to_dict.side_effect = lambda: self.header
        self.bam.__enter__.return_value = self.bam
        self.pysam = MagicMock(); self.pysam.AlignmentFile.return_value = self.bam

    def validate(self):
        with patch.dict(sys.modules, {'pysam': self.pysam}):
            return braker.validate_inputs(self.spec, self.genome, self.protein,
                                         ['reads.bam'], [self.index], self.identity)

    def test_checks_correct_dictionary_readgroup_and_explicit_csi(self):
        self.assertEqual(self.validate()['status'], 'ready')
        self.pysam.AlignmentFile.assert_called_once_with('reads.bam', 'rb', index_filename=str(self.index))
        self.bam.check_index.assert_called_once()

    def test_changed_bytes_dictionary_and_wrong_readgroup_block(self):
        self.protein.write_text('>other\nMFK\n')
        with self.assertRaisesRegex(ValueError, 'bytes differ'): self.validate()
        self.protein.write_text('>p\nMFK\n')
        self.bam.lengths = (9,)
        with self.assertRaisesRegex(ValueError, 'dictionary'): self.validate()
        self.bam.lengths = (8,); self.header['RG'][0]['SM'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'wrong readset'): self.validate()


if __name__ == '__main__':
    unittest.main()
