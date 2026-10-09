"""Reference and repeat integrity contracts; executed by the cluster acceptance runner."""
import csv
import gzip
import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'py_scripts'))
import reference_repeats as rr
from report_reference_repeats import summarize


class ReferenceRepeatsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_embedded_scripts_escape_literal_dollars_for_groovy(self):
        # Python/shell quoting does not suppress interpolation in a Groovy GString.
        pending, seen = [ROOT / 'reference_repeats.nf'], set()
        while pending:
            path = pending.pop().resolve()
            if path in seen:
                continue
            seen.add(path)
            text = path.read_text()
            for relative in re.findall(r"\bfrom\s+['\"]([^'\"]+)['\"]", text):
                pending.append((path.parent / relative).with_suffix('.nf'))
            for block in re.finditer('"""(.*?)"""', text, re.S):
                script = block[1]
                for position, character in enumerate(script):
                    if character != '$':
                        continue
                    slashes, previous = 0, position - 1
                    while previous >= 0 and script[previous] == '\\':
                        slashes += 1
                        previous -= 1
                    if slashes % 2 == 0 and script[position + 1:position + 2] != '{':
                        line = text[:block.start(1) + position].count('\n') + 1
                        self.fail(f'{path.name}:{line}: escape literal dollars or use braced Groovy interpolation')

    def write(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def prepare(self, text='>long__id_100-500 description\n' + 'ACGT' * 40 + '\n>small\nACGT\n', minimum=100):
        source = self.write('source.fa', text)
        output = self.root / 'reference'
        state = rr.reference(source, 'A', minimum, output)
        return state, output / 'A.filtered.fasta', output / 'A.reference.json'

    def hit(self, start, end, linkage='1', name='long__id_100-500', family='DNA/hAT'):
        return f'500 0.0 0.0 0.0 {name} {start} {end} (0) + repeat {family} 1 30 (0) {linkage}\n'

    def test_reference_threshold_index_and_reversible_headers(self):
        state, fasta, _ = self.prepare('>unsafe|name full original header\n' + 'acgt' * 40 + '\n>small\nACGT\n')
        self.assertEqual(state['status'], 'ready')
        self.assertEqual(state['stats']['excluded_bases'], 4)
        self.assertEqual(list(rr.fasta(fasta))[0][0], 'unsafe_name')
        mapping = list(csv.DictReader((fasta.parent / 'A.header_map.tsv').open(), delimiter='\t'))
        self.assertEqual(mapping[0]['original_header'], 'unsafe|name full original header')
        fields = (fasta.parent / 'A.filtered.fasta.fai').read_text().split()
        with fasta.open('rb') as handle:
            handle.seek(int(fields[2]))
            self.assertEqual(handle.read(int(fields[3])), ('ACGT' * 15).encode())
        self.assertEqual(fields[1:], ['160', '13', '60', '61'])

    def test_zero_threshold_and_gzip(self):
        source = self.root / 'genome.gz'
        with gzip.open(source, 'wt') as handle:
            handle.write('>a\nacgt\n>b\nN\n')
        state = rr.reference(source, 'G', 0, self.root / 'gz')
        self.assertEqual(state['stats']['retained_contigs'], 2)
        self.assertEqual(state['stats']['excluded_contigs'], 0)

    def test_threshold_is_inclusive(self):
        state, _, _ = self.prepare(minimum=160)
        self.assertEqual(state['status'], 'ready')
        self.assertEqual(state['stats']['retained_contigs'], 1)

    def test_collisions_even_excluded_and_invalid_sequences(self):
        for index, text in enumerate(('>a|b\nACGT\n>a_b\n' + 'A' * 160 + '\n',
                                      '>a\nAC-GT\n', '>a\nACZ\n', '>a\n', 'ACGT\n', '>a\nAC GT\n')):
            source = self.write(f'bad{index}.fa', text)
            folder = self.root / f'bad{index}'
            self.assertEqual(rr.reference(source, 'A', 100, folder)['status'], 'failed')
            self.assertFalse((folder / 'A.filtered.fasta').exists())

    def test_record_safe_curation(self):
        source = self.write('families.fa', '>simple#Simple_repeat\nAAAA\nCCCC\n>keep#DNA/hAT\nAC-G\nTAC.G\n'
                            '>duplicate#Unknown\nACGTACG\n>short#Unknown\nAC\n')
        output, report = self.root / 'curated.fa', self.root / 'curation.tsv'
        rr.curate(source, output, 4, report)
        self.assertEqual(list(rr.fasta(output)), [('keep#DNA/hAT', 'ACGTACG')])
        self.assertEqual([row['decision'] for row in csv.DictReader(report.open(), delimiter='\t')],
                         ['excluded_class', 'retained', 'duplicate_sequence', 'short'])

    def test_whole_contigs_no_empty_batches_and_coordinate_helpers(self):
        state, fasta, _ = self.prepare(minimum=0)
        folder = self.root / 'batches'
        rr.batches(fasta, 'A', 8, folder)
        accounting = json.loads((folder / 'A.batches.json').read_text())
        self.assertEqual(accounting['actual_batches'], 2)
        self.assertEqual(sum(accounting['batch_bases']), state['stats']['retained_bases'])
        records = [row for path in sorted((folder / 'batches').glob('*.fa')) for row in rr.fasta(path)]
        self.assertCountEqual(records, list(rr.fasta(fasta)))
        for index, path in enumerate(sorted((folder / 'batches').glob('*.fa'))):
            mapping = self.root / f'map{index}.tsv'
            subprocess.run([sys.executable, str(ROOT / 'py_scripts/build_coord_map.py'), str(path),
                            str(folder / 'A_seq_lengths.tsv'), str(mapping)], check=True, capture_output=True)
            for row in csv.reader(mapping.open(), delimiter='\t'):
                self.assertEqual(row[2], '1')
                self.assertEqual(row[3], row[4])
        restored = self.write('restored.fa', '>long__id_100-500__1__160\n' + 'ACGT' * 40 + '\n')
        subprocess.run([sys.executable, str(ROOT / 'py_scripts/verify_coord_map.py'), str(restored)],
                       check=True, capture_output=True)

    def test_merge_restores_order_and_long_header_without_folding(self):
        name = 'contig__' + 'x' * 75
        lengths = self.write('lengths.tsv', f'{name}\t4\nb\t4\n')
        first = self.write('one.fa', '>b__1__4\naCGT\n')
        second = self.write('two.fa', f'>{name}__1__4\nACgt\n')
        output = self.root / 'merged.fa'
        rr.merge_fasta([first, second], lengths, output)
        self.assertEqual(list(rr.fasta(output)), [(name, 'ACgt'), ('b', 'aCGT')])
        self.assertEqual(output.read_text().splitlines()[0], '>' + name)
        for files in ([first], [first, second, second]):
            with self.assertRaises(ValueError):
                rr.merge_fasta(files, lengths, output)
        second.write_text(f'>{name}__2__4\nCGT\n')
        with self.assertRaises(ValueError):
            rr.merge_fasta([first, second], lengths, output)

    def test_linked_fragments_and_overlapping_coverage(self):
        state, fasta, reference = self.prepare()
        lengths = self.write('lengths.tsv', 'long__id_100-500\t160\n')
        one = self.write('one.out', self.hit(1, 30) + self.hit(20, 50))
        two = self.write('two.out', self.hit(60, 80))
        merged = self.root / 'merged.out'
        rr.merge_out([one, two], lengths, merged)
        self.assertEqual([row[14] for row in rr.out_rows(merged)], ['1', '1', '2'])
        masked = self.write('masked.fa', '>long__id_100-500\n' + ('ACGT' * 40).lower() + '\n')
        result = rr.validate_mask(reference, masked, merged, fasta, 'A', self.root / 'result.json')
        self.assertEqual(result['status'], 'ready')
        self.assertEqual(result['repeat_interval_union_bases'], 71)
        self.assertEqual(result['repeat_alignment_records'], 3)
        self.assertEqual(result['linked_repeat_elements'], 2)
        self.assertEqual(result['lowercase_acgt_bases'], 160)
        self.assertEqual(rr.union_length([(1, 10), (5, 20), (21, 21), (30, 31)]), 23)

    def test_mask_identity_and_coordinate_failures(self):
        state, fasta, reference = self.prepare()
        empty = self.write('empty.out', 'header\n')
        self.assertEqual(rr.validate_mask(reference, fasta, empty, fasta, 'A', self.root / 'ok.json')['status'], 'ready')
        wrong = self.write('wrong.fa', '>long__id_100-500\n' + 'A' * 160 + '\n')
        self.assertEqual(rr.validate_mask(reference, wrong, empty, fasta, 'A', self.root / 'bad.json')['status'], 'failed')
        bad_out = self.write('bad.out', self.hit(1, 161))
        self.assertEqual(rr.validate_mask(reference, fasta, bad_out, fasta, 'A', self.root / 'bounds.json')['status'], 'failed')

    def test_library_failure_and_unsplit_comparison(self):
        state, fasta, reference = self.prepare()
        for index, text in enumerate(('', '>x\nACGT\n>x\nTGCA\n')):
            path = self.write(f'lib{index}.fa', text)
            self.assertEqual(rr.validate_library(path, 'A', self.root / 'lib.json')['status'], 'failed')
        one = self.write('one.out', self.hit(1, 30, '1'))
        two = self.write('two.out', self.hit(1, 30, '999'))
        self.assertEqual(rr.compare_masks(fasta, one, [fasta, two], 'A', self.root / 'benchmark.json')['status'], 'ready')
        two.write_text(self.hit(2, 30))
        self.assertEqual(rr.compare_masks(fasta, one, [fasta, two], 'A', self.root / 'benchmark.json')['status'], 'failed')

    def test_report_keeps_independent_ready_assembly(self):
        report = summarize({'assemblies': [{'assembly_id': 'A'}, {'assembly_id': 'B'}]},
                           [{'assembly_id': 'A', 'status': 'ready'}, {'assembly_id': 'B', 'status': 'failed', 'reason': 'bad DNA'}],
                           [{'assembly_id': 'A', 'status': 'ready'}], [{'assembly_id': 'A', 'status': 'ready'}], self.root)
        self.assertEqual(report['status'], 'incomplete')
        self.assertEqual([row['status'] for row in report['assemblies']], ['ready', 'blocked'])


if __name__ == '__main__':
    unittest.main()
