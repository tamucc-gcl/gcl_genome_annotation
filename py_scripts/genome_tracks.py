#!/usr/bin/env python3
"""Validate independent genome tracks and summarize them without changing gene models."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from urllib.parse import quote

from reference_repeats import fasta, sha, dump

TRACKS = ('trnascan', 'infernal', 'telomere', 'centromere')


def dictionary(path):
    data = json.loads(Path(path).read_text())
    return data, {r['id']: r['length'] for r in data['contigs']}


def centromeres(monomers, hors, aid, output):
    """centroAnno 1.0.2 writers use zero-based INCLUSIVE start/end coordinates."""
    lines = ['##gff-version 3\n']
    for source, kind, width in ((monomers, 'satellite_DNA', 6), (hors, 'repeat_region', 7)):
        with Path(source).open() as handle:
            for cols in csv.reader(handle):
                if not cols or cols[0].startswith('#') or cols[0].strip().lower() == 'sequence name':
                    continue
                if len(cols) != width:
                    raise ValueError('Unexpected centroAnno CSV schema: ' + str(cols))
                seq, template = cols[:2]
                start, end, size = int(cols[2]), int(cols[3]), int(cols[5])
                identity = float(cols[4])
                if start < 0 or end < start or end - start + 1 != size or not 0 <= identity <= 1:
                    raise ValueError('Invalid centroAnno coordinates/length/identity: ' + str(cols))
                attrs = f'ID={aid}_centro{len(lines)};Name={quote(template, safe="._-")};length={size};identity={identity}'
                if width == 7:
                    attrs += ';confidence=' + quote(cols[6], safe='._-')
                lines.append(f'{seq}\tcentroAnno\t{kind}\t{start+1}\t{end+1}\t{identity*1000:.1f}\t.\t.\t{attrs}\n')
    Path(output).write_text(''.join(lines))


def telomeres(tsv, genome, aid, motif, window, minimum, terminal, output):
    """tidk TSV uses window END coordinates, including the clipped final window."""
    if window < 1 or minimum < 1 or terminal < 0:
        raise ValueError('Invalid telomere window/count/end distance')
    lengths = {h.split()[0]: len(s) for h, s in fasta(genome)}
    grouped = defaultdict(list)
    with Path(tsv).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            seq = row['id']
            end = int(row['window'])
            if seq not in lengths:
                raise ValueError('Unknown tidk contig: ' + seq)
            start = len(grouped[seq]) * window
            if end != min(start + window, lengths[seq]) or start >= end:
                raise ValueError('Unexpected tidk window end/order: ' + str(row))
            fwd, rev = int(row['forward_repeat_number']), int(row['reverse_repeat_number'])
            if min(fwd, rev) < 0 or row['telomeric_repeat'].upper() != motif.upper():
                raise ValueError('Invalid tidk count/motif')
            grouped[seq].append((start, end, fwd, rev))
    if set(grouped) != set(lengths) or any(rows[-1][1] != lengths[s] for s, rows in grouped.items()):
        raise ValueError('Incomplete tidk windows')
    lines = ['##gff-version 3\n']
    for seq in sorted(grouped):
        regions = []
        for start, end, fwd, rev in grouped[seq]:
            if fwd + rev < minimum:
                continue
            if regions and regions[-1][1] == start:
                regions[-1][1] = end
                regions[-1][2] += fwd
                regions[-1][3] += rev
            else:
                regions.append([start, end, fwd, rev])
        for start, end, fwd, rev in regions:
            left, right = start <= terminal, end >= lengths[seq] - terminal
            position = 'spanning' if left and right else '5_prime' if left else '3_prime' if right else 'interstitial'
            attrs = f'ID={aid}_telo{len(lines)};Name=telomeric_repeat;motif={motif};n_repeats={fwd+rev};forward_count={fwd};reverse_count={rev};position={position}'
            lines.append(f'{seq}\ttidk\ttelomeric_repeat\t{start+1}\t{end}\t{fwd+rev}\t{"+" if fwd >= rev else "-"}\t.\t{attrs}\n')
    Path(output).write_text(''.join(lines))


def normalize(aid, track, source, identity, output):
    data, lengths = dictionary(identity)
    rows, ids, parents = [], set(), []
    with Path(source).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) != 9:
                raise ValueError('Expected nine GFF columns')
            start, end = int(r[3]), int(r[4])
            if r[0] not in lengths or not 1 <= start <= end <= lengths[r[0]] or r[6] not in ('+', '-', '.', '?'):
                raise ValueError('Track coordinates/strand outside reference: ' + line.strip())
            attrs = dict(x.strip().split('=', 1) for x in r[8].strip(';').split(';') if '=' in x)
            # Infernal's existing converter uses the Rfam family as column 3.
            if track == 'infernal':
                attrs['Name'] = quote(r[2], safe='._-')
                r[2] = 'ncRNA'
            old = attrs.get('ID', f'feature{len(rows)+1}')
            attrs['ID'] = f'{aid}.{track}.{old}'
            if attrs['ID'] in ids:
                raise ValueError('Duplicate feature ID: ' + old)
            ids.add(attrs['ID'])
            if 'Parent' in attrs:
                attrs['Parent'] = ','.join(f'{aid}.{track}.{p}' for p in attrs['Parent'].split(','))
                parents.extend(attrs['Parent'].split(','))
            r[8] = ';'.join(f'{k}={v}' for k, v in attrs.items())
            rows.append(r)
    if set(parents) - ids:
        raise ValueError('Missing track feature parents')
    rows.sort(key=lambda r: (r[0], int(r[3]), int(r[4]), r[2], r[8]))
    Path(output).write_text('##gff-version 3\n' + ''.join('\t'.join(r) + '\n' for r in rows))
    state = dict(status='ready', assembly_id=aid, track=track, features=len(rows),
                 feature_types=dict(Counter(r[2] for r in rows)),
                 reference_identity=data['reference_identity'], source_sha256=sha(source), gff3_sha256=sha(output))
    dump(track + '.status.json', state)
    return state



def main():
    parser = argparse.ArgumentParser(description=__doc__); subs = parser.add_subparsers(dest='action', required=True)
    p = subs.add_parser('centromeres')
    for name in ('monomers', 'hors', 'aid', 'output'): p.add_argument('--'+name, required=True)
    p = subs.add_parser('telomeres')
    for name in ('tsv', 'genome', 'aid', 'motif', 'output'): p.add_argument('--'+name, required=True)
    for name in ('window', 'minimum', 'terminal'): p.add_argument('--'+name, required=True, type=int)
    p = subs.add_parser('normalize')
    for name in ('aid', 'source', 'identity', 'output'): p.add_argument('--'+name, required=True)
    p.add_argument('--track', required=True, choices=TRACKS)
    a = vars(parser.parse_args()); action = a.pop('action')
    globals()[action](**a)


if __name__ == '__main__': main()
