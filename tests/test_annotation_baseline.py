import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from annotation_baseline import prepare, report, sha


def fixture(folder, aid='A'):
    folder = Path(folder)
    gff, proteins = folder / f'{aid}.gff3', folder / f'{aid}.fa'
    # Features deliberately out of parent order; shared CDS has two parents.
    gff.write_text('##gff-version 3\n' + '\n'.join([
        's\tBRAKER\tCDS\t1\t30\t.\t+\t0\tParent=t2,t1',
        's\tBRAKER\tmRNA\t1\t30\t.\t+\t.\tID=t2;Parent=g1',
        's\tBRAKER\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1',
        's\tBRAKER\tgene\t1\t30\t.\t+\t.\tID=g1',
        's\tBRAKER\tgene\t40\t80\t.\t+\t.\tID=g2',
        's\tBRAKER\tmRNA\t40\t80\t.\t+\t.\tID=t3;Parent=g2',
        's\tBRAKER\tCDS\t40\t80\t.\t+\t0\tParent=t3']) + '\n')
    proteins.write_text('>t2\nMAAAA*\n>t1\nMAAAA\n>t3\nMXX*AA\n')
    status = dict(assembly_id=aid, status='ready', genes=2, coding_transcripts=3,
                  gff3_sha256=sha(gff), proteins_sha256=sha(proteins))
    return gff, proteins, status


class BaselineContracts(unittest.TestCase):
    def test_mismatched_genome_blocks_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff, proteins, status = fixture(folder)
            genome = folder / 'genome.fa'; genome.write_text('>s\n' + 'A' * 100 + '\n')
            with self.assertRaisesRegex(ValueError, 'Genome differs'):
                prepare('A', {}, status | {'reference_identity': 'different'}, gff, proteins, folder / 'out', genome)

    def test_representatives_explicit_parents_stable_ties_no_model_filtering(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff, proteins, status = fixture(folder)
            before = (gff.read_bytes(), proteins.read_bytes())
            row = prepare('A', {'run_busco': False, 'run_omark': False}, status, gff, proteins, folder / 'out')
            self.assertEqual(row['representative_proteins'], 2)
            self.assertEqual(row['coding_transcripts'], 3)
            self.assertEqual((folder / 'out/A.fa').read_text(), '>t1\nMAAAA\n>t3\nMXX*AA\n')
            self.assertEqual(before, (gff.read_bytes(), proteins.read_bytes()))

    def test_changed_predictions_wrong_assembly_unvalidated_outputs_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff, proteins, status = fixture(folder)
            for update in ({'assembly_id': 'B'}, {'status': 'software_checked'}, {'proteins_sha256': 'wrong'}, {'genes': 9}):
                with self.subTest(update=update), self.assertRaises(ValueError):
                    prepare('A', {}, status | update, gff, proteins, folder / 'out')

    def test_requested_metrics_require_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff, proteins, status = fixture(tmp)
            for settings in ({'run_busco': True}, {'run_omark': True}):
                with self.subTest(settings=settings), self.assertRaises(ValueError):
                    prepare('A', settings, status, gff, proteins, Path(tmp) / 'out')

    def test_disabled_metrics_report_without_placeholder_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff, proteins, status = fixture(folder)
            prepare('A', {}, status, gff, proteins, folder / 'selection')
            stats = folder / 'agat.txt'; stats.write_text('actual AGAT stats\n')
            row = report([folder / 'selection/selection.json'], [{'aid': 'A', 'path': stats}], [], [], folder / 'report')
            self.assertEqual(row['assemblies'][0]['analyses'], dict(agat='ready', busco='not_requested', omark='not_requested'))
            self.assertNotIn('BUSCO summary', (folder / 'report/README.md').read_text())

    def test_missing_requested_or_empty_metrics_and_duplicate_assemblies_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            gff, proteins, status = fixture(folder)
            prepare('A', {'run_busco': True, 'busco_lineage': 'actinopterygii_odb10'}, status, gff, proteins, folder / 'selection')
            stats = folder / 'agat.txt'; stats.write_text('AGAT stats\n')
            selections = [folder / 'selection/selection.json']
            with self.assertRaisesRegex(ValueError, 'Missing requested busco'):
                report(selections, [{'aid': 'A', 'path': stats}], [], [], folder / 'report')
            with self.assertRaises(ValueError):
                report(selections, [{'aid': 'A', 'path': stats}], [{'aid': 'B', 'path': stats}], [], folder / 'report')
            with self.assertRaisesRegex(ValueError, 'Duplicate baseline'):
                report(selections * 2, [{'aid': 'A', 'path': stats}], [{'aid': 'A', 'path': stats}], [], folder / 'report')
            stats.write_text('')
            with self.assertRaisesRegex(ValueError, 'Empty/missing'):
                report(selections, [{'aid': 'A', 'path': stats}], [], [], folder / 'report')


if __name__ == '__main__':
    unittest.main()
