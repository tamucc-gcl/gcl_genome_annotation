#!/usr/bin/env python3
"""Reference validation and record-safe repeat preparation; Nextflow owns file staging."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

DNA = set('ACGTRYSWKMBDHVN')


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fasta(path, allow_gaps=False):
    with Path(path).open('rb') as handle:
        zipped = handle.read(2) == b'\x1f\x8b'
    opener = gzip.open if zipped else open
    header, lines = None, []
    with opener(path, 'rt', encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    if not lines:
                        raise ValueError(f'Empty FASTA sequence: {header}')
                    yield header, ''.join(lines)
                header, lines = line[1:].strip(), []
                if not header:
                    raise ValueError(f'Empty FASTA header at line {number}')
            else:
                if header is None:
                    raise ValueError(f'Sequence before FASTA header at line {number}')
                if any(char.isspace() for char in line):
                    raise ValueError(f'Whitespace inside sequence at line {number}')
                if not set(line.upper()) <= (DNA | set('-.') if allow_gaps else DNA):
                    raise ValueError(f'Invalid nucleotide/gap at line {number}; coordinates cannot be repaired silently')
                lines.append(line)
    if header is not None:
        if not lines:
            raise ValueError(f'Empty FASTA sequence: {header}')
        yield header, ''.join(lines)


def write_record(handle, header, sequence):
    handle.write('>' + header + '\n')
    for start in range(0, len(sequence), 60):
        handle.write(sequence[start:start + 60] + '\n')


def reference(source, aid, minimum, outdir):
    """Normalize IDs and bases, retain >= minimum, emit dictionary and valid faidx offsets."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / aid
    status = {'schema_version': 1, 'stage': 'reference', 'assembly_id': aid,
              'status': 'failed', 'min_contig_length': minimum, 'implementation_sha256': sha(__file__)}
    try:
        if minimum < 0:
            raise ValueError('Minimum contig length must be nonnegative')
        rows, seen, kept, excluded = [], set(), [], []
        with Path(str(prefix) + '.filtered.fasta').open('wb') as output:
            for header, sequence in fasta(source):
                original = header.split()[0]
                clean = re.sub('[^A-Za-z0-9._-]', '_', original)
                if clean in seen:
                    raise ValueError(f'Duplicate/colliding normalized contig ID: {clean}')
                seen.add(clean)
                sequence = sequence.upper()
                retain = len(sequence) >= minimum
                rows.append([clean, original, header, 'retained' if retain else 'removed', len(sequence)])
                record = {'id': clean, 'length': len(sequence),
                          'uppercase_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest()}
                if not retain:
                    excluded.append(record)
                    continue
                output.write(('>' + clean + '\n').encode('ascii'))
                offset = output.tell()
                for start in range(0, len(sequence), 60):
                    output.write((sequence[start:start + 60] + '\n').encode('ascii'))
                record['fai'] = [clean, len(sequence), offset, min(60, len(sequence)), min(60, len(sequence)) + 1]
                kept.append(record)
        if not kept:
            raise ValueError('No contigs survived the configured minimum length')
        with Path(str(prefix) + '.header_map.tsv').open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['clean_id', 'original_id', 'original_header', 'status', 'length'])
            writer.writerows(rows)
        with Path(str(prefix) + '.excluded_contigs.tsv').open('w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['contig', 'length', 'reason'])
            writer.writerows([row['id'], row['length'], 'below_min_contig_length'] for row in excluded)
        with Path(str(prefix) + '.filtered.fasta.fai').open('w', encoding='utf-8') as handle:
            for row in kept:
                handle.write('\t'.join(map(str, row.pop('fai'))) + '\n')
        dictionary = {'schema_version': 1, 'assembly_id': aid, 'contigs': kept}
        identity = hashlib.sha256(json.dumps(kept, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        dictionary['reference_identity'] = identity
        dump(str(prefix) + '.contigs.json', dictionary)
        stats = {'input_contigs': len(rows), 'input_bases': sum(row[4] for row in rows),
                 'retained_contigs': len(kept), 'retained_bases': sum(row['length'] for row in kept),
                 'excluded_contigs': len(excluded), 'excluded_bases': sum(row['length'] for row in excluded)}
        with Path(str(prefix) + '.assembly_stats.tsv').open('w', encoding='utf-8') as handle:
            handle.write('metric\tvalue\n' + ''.join(f'{key}\t{value}\n' for key, value in stats.items()))
        status.update(status='ready', source_sha256=sha(source), reference_identity=identity,
                      normalized_sha256=sha(str(prefix) + '.filtered.fasta'), stats=stats,
                      excluded_lengths=[row['length'] for row in excluded], contigs=kept)
    except (ValueError, OSError, EOFError, UnicodeError) as exc:
        status['reason'] = str(exc)
        # Partial scientific outputs must never be advertised as accepted references.
        for suffix in ('.filtered.fasta', '.filtered.fasta.fai', '.contigs.json'):
            Path(str(prefix) + suffix).unlink(missing_ok=True)
    dump(str(prefix) + '.reference.json', status)
    return status


def curate(source, output, minimum, report):
    seen, rows = set(), []
    with Path(output).open('w', encoding='utf-8') as handle:
        for header, sequence in fasta(source, allow_gaps=True):
            sequence = sequence.replace('-', '').replace('.', '')
            classification = header.split()[0].partition('#')[2].split('/')[0]
            reason = ('short' if len(sequence) < minimum else
                      'excluded_class' if classification in ('Simple_repeat', 'Low_complexity', 'Satellite') else
                      'duplicate_sequence' if sequence.upper() in seen else 'retained')
            rows.append([header.split()[0], len(sequence), reason])
            if reason == 'retained':
                seen.add(sequence.upper())
                write_record(handle, header, sequence)
    with Path(report).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['repeat_id', 'length', 'decision'])
        writer.writerows(rows)


def batches(source, aid, count, outdir):
    """Balance whole contigs across at most count batches; never tile a contig."""
    if count < 1:
        raise ValueError('n_rm_chunks must be positive')
    lengths = []
    for header, sequence in fasta(source):
        if header != header.split()[0] or not re.fullmatch('[A-Za-z0-9._-]+', header):
            raise ValueError('Masking batches require normalized bare FASTA IDs')
        lengths.append((header, len(sequence)))
    if not lengths or len({row[0] for row in lengths}) != len(lengths):
        raise ValueError('Empty reference or duplicate contig ID')
    n = min(count, len(lengths))
    loads, assignments = [0] * n, {}
    for name, size in sorted(lengths, key=lambda row: (-row[1], row[0])):
        index = min(range(n), key=lambda index: (loads[index], index))
        assignments[name] = index
        loads[index] += size
    outdir = Path(outdir)
    (outdir / 'batches').mkdir(parents=True)
    (outdir / 'beds').mkdir()
    handles = [(outdir / f'batches/batch_{i:04d}.fa').open('w', encoding='utf-8') for i in range(n)]
    beds = [(outdir / f'beds/batch_{i:04d}.bed').open('w', encoding='utf-8') for i in range(n)]
    try:
        for header, sequence in fasta(source):
            index = assignments[header]
            write_record(handles[index], header, sequence)
            beds[index].write(f'{header}\t0\t{len(sequence)}\n')
    finally:
        for handle in handles + beds:
            handle.close()
    (outdir / f'{aid}_seq_lengths.tsv').write_text(''.join(f'{name}\t{size}\n' for name, size in lengths))
    dump(outdir / f'{aid}.batches.json', {'schema_version': 1, 'policy': 'whole_contigs', 'requested_batches': count,
                                       'actual_batches': n, 'batch_bases': loads, 'contigs': len(lengths)})


def merge_fasta(files, lengths, output):
    expected = dict((name, int(size)) for name, size in csv.reader(Path(lengths).open(), delimiter='\t'))
    sequences = {}
    for path in files:
        for header, sequence in fasta(path):
            match = re.fullmatch(r'(.+)__(\d+)__(\d+)', header)
            if not match:
                raise ValueError(f'Missing restored FASTA coordinates: {header}')
            name, start, end = match[1], int(match[2]), int(match[3])
            if name not in expected or start != 1 or end != expected[name] or len(sequence) != expected[name]:
                raise ValueError(f'Incomplete/incorrect whole-contig reconstruction: {header}')
            if name in sequences:
                raise ValueError(f'Duplicate reconstructed contig: {name}')
            sequences[name] = sequence
    if set(sequences) != set(expected):
        raise ValueError('Masked contig dictionary does not match the input reference')
    with Path(output).open('w', encoding='utf-8') as handle:
        for name in expected:
            write_record(handle, name, sequences[name])


def out_rows(path):
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        if len(fields) < 15:
            raise ValueError(f'Malformed RepeatMasker row {path}:{number}')
        yield fields


def merge_out(files, lengths, output):
    expected = dict((name, int(size)) for name, size in csv.reader(Path(lengths).open(), delimiter='\t'))
    linkage, counter = {}, 0
    with Path(output).open('w', encoding='utf-8') as handle:
        handle.write('   SW  perc perc perc  query       position in query    matching repeat\n'
                     'score  div. del. ins.  sequence    begin end (left)     repeat class/family\n\n')
        for path in sorted(map(Path, files)):
            for fields in out_rows(path):
                name, start, end = fields[4], int(fields[5]), int(fields[6])
                if name not in expected or not 1 <= start <= end <= expected[name]:
                    raise ValueError(f'Out-of-bounds repeat: {name}:{start}-{end}')
                key = (str(path), fields[14])
                if key not in linkage:
                    counter += 1
                    linkage[key] = counter
                fields[14] = str(linkage[key])
                fields[7] = f'({expected[name] - end})'
                handle.write(' '.join(fields) + '\n')


def union_length(intervals):
    total, left, right = 0, None, None
    for start, end in sorted(intervals):
        if left is None:
            left, right = start, end
        elif start <= right + 1:
            right = max(right, end)
        else:
            total += right - left + 1
            left, right = start, end
    return total + (right - left + 1 if left is not None else 0)


def validate_library(source, aid, output):
    result = {'schema_version': 1, 'assembly_id': aid, 'status': 'failed'}
    try:
        seen, bases = set(), 0
        for header, sequence in fasta(source):
            name = header.split()[0]
            if name in seen:
                raise ValueError(f'Duplicate repeat-library ID: {name}')
            seen.add(name)
            bases += len(sequence)
        if not seen:
            raise ValueError('Combined repeat library is empty')
        result.update(status='ready', library_sha256=sha(source), families=len(seen), bases=bases)
    except (ValueError, OSError, EOFError, UnicodeError) as exc:
        result['reason'] = str(exc)
    dump(output, result)
    return result


def compare_masks(masked, annotations, files, aid, output):
    try:
        if list(fasta(masked)) != list(fasta(files[0])):
            raise ValueError('Masked sequences differ from the unsplit reference run')
        # Ignore globally reassigned linkage IDs and output row ordering.
        scientific_hits = lambda path: sorted(tuple(row[4:7] + row[8:14]) for row in out_rows(path))
        if scientific_hits(annotations) != scientific_hits(files[1]):
            raise ValueError('Repeat positions/classes/targets differ from the unsplit reference run')
        result = {'assembly_id': aid, 'status': 'ready', 'comparison': 'exact masked bases and scientific hit coordinates'}
    except (ValueError, OSError, EOFError, UnicodeError) as exc:
        result = {'assembly_id': aid, 'status': 'failed', 'reason': str(exc)}
    dump(output, result)
    return result


def validate_mask(reference_json, masked, annotations, library, aid, output):
    result = {'schema_version': 1, 'assembly_id': aid, 'stage': 'reference_repeats', 'status': 'failed'}
    try:
        reference_data = json.loads(Path(reference_json).read_text())
        expected = reference_data['contigs']
        observed, seen, lowercase = [], set(), 0
        for name, sequence in fasta(masked):
            if name in seen:
                raise ValueError('Duplicate masked contig')
            seen.add(name)
            lowercase += sum(char in 'acgt' for char in sequence)
            observed.append({'id': name, 'length': len(sequence),
                             'uppercase_sha256': hashlib.sha256(sequence.upper().encode('ascii')).hexdigest()})
        if observed != expected:
            raise ValueError('Masking changed reference contig order/IDs, lengths, or uppercase bases')
        intervals, classes, elements, hits = defaultdict(list), defaultdict(lambda: defaultdict(list)), set(), 0
        lengths = {row['id']: row['length'] for row in expected}
        for row in out_rows(annotations):
            name, start, end = row[4], int(row[5]), int(row[6])
            if name not in lengths or not 1 <= start <= end <= lengths[name]:
                raise ValueError('Repeat annotation outside normalized reference')
            intervals[name].append((start, end))
            classes[row[10]][name].append((start, end))
            elements.add(row[14])
            hits += 1
        coverage = sum(union_length(value) for value in intervals.values())
        class_bp = {key: sum(union_length(value) for value in groups.values()) for key, groups in classes.items()}
        result.update(status='ready', reference_identity=reference_data['reference_identity'],
                      masked_sha256=sha(masked), repeat_out_sha256=sha(annotations), library_sha256=sha(library),
                      contigs=len(expected), reference_bases=sum(lengths.values()), lowercase_acgt_bases=lowercase,
                      repeat_interval_union_bases=coverage, repeat_alignment_records=hits,
                      linked_repeat_elements=len(elements), class_union_bases=class_bp,
                      statistics_policy='Inclusive interval unions; class unions overlap and must not be summed')
    except (ValueError, OSError, KeyError, EOFError, UnicodeError) as exc:
        result['reason'] = str(exc)
    dump(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['reference', 'curate', 'batches', 'merge-fasta', 'merge-out', 'library', 'compare', 'validate'])
    for name in ('source', 'assembly', 'outdir', 'output', 'report', 'lengths', 'reference-json', 'masked', 'annotations', 'library'):
        parser.add_argument('--' + name)
    parser.add_argument('--minimum', type=int, default=1000)
    parser.add_argument('--count', type=int, default=8)
    parser.add_argument('--files', nargs='+')
    args = parser.parse_args()
    if args.action == 'reference':
        reference(args.source, args.assembly, args.minimum, args.outdir)
    elif args.action == 'curate':
        curate(args.source, args.output, args.minimum, args.report)
    elif args.action == 'batches':
        batches(args.source, args.assembly, args.count, args.outdir)
    elif args.action == 'merge-fasta':
        merge_fasta(args.files, args.lengths, args.output)
    elif args.action == 'merge-out':
        merge_out(args.files, args.lengths, args.output)
    elif args.action == 'library':
        validate_library(args.source, args.assembly, args.output)
    elif args.action == 'compare':
        compare_masks(args.masked, args.annotations, args.files, args.assembly, args.output)
    else:
        validate_mask(args.reference_json, args.masked, args.annotations, args.library, args.assembly, args.output)


if __name__ == '__main__':
    main()
