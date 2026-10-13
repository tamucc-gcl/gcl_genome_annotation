#!/usr/bin/env python3
"""Validate the unfiltered model handoff and select one protein per coding gene.

No model filtering or changes to the accepted prediction files are performed.
"""
import argparse
import base64
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote

from prepare_protein_evidence import records, checksum as sha
from reference_repeats import fasta as genome_records


def prepare(aid, settings, prediction, gff, proteins, output, genome=None):
    if prediction.get('assembly_id') != aid or prediction.get('status') != 'ready':
        raise ValueError('Baseline requires validated predictions for the requested assembly')
    if sha(gff) != prediction['gff3_sha256'] or sha(proteins) != prediction['proteins_sha256']:
        raise ValueError('Prediction files changed since validation')
    if genome is not None:
        contigs = [{'id': header.split()[0], 'length': len(sequence),
                    'uppercase_sha256': hashlib.sha256(sequence.upper().encode('ascii')).hexdigest()}
                   for header, sequence in genome_records(genome)]
        identity = hashlib.sha256(json.dumps(contigs, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        if identity != prediction.get('reference_identity'):
            raise ValueError('Genome differs from the accepted prediction reference')
    if settings.get('run_busco') and not settings.get('busco_lineage'):
        raise ValueError('Enabled BUSCO requires a resolved lineage')
    if settings.get('run_omark') and not settings.get('taxid'):
        raise ValueError('Enabled OMArk requires a resolved taxid')
    genes, parents, coding = set(), {}, set()
    with Path(gff).open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('Malformed prediction GFF3')
            attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            identifier = unquote(attrs.get('ID', ''))
            links = [unquote(x) for x in attrs.get('Parent', '').split(',') if x]
            if fields[2] == 'gene':
                if not identifier or identifier in genes:
                    raise ValueError('Missing/duplicate gene ID')
                genes.add(identifier)
            elif fields[2] in ('mRNA', 'transcript'):
                if not identifier or identifier in parents or len(links) != 1:
                    raise ValueError('Transcript requires one explicit gene parent')
                parents[identifier] = links[0]
            elif fields[2] == 'CDS':
                coding.update(links)
    if not coding or not coding <= parents.keys() or not set(parents.values()) <= genes:
        raise ValueError('Unresolved coding transcript/gene relationships')
    sequences = {}
    with Path(proteins).open(encoding='utf-8') as handle:
        for header, sequence in records(handle):
            tid = header.split()[0]
            if tid in sequences or not sequence:
                raise ValueError('Duplicate/empty protein sequence')
            sequences[tid] = sequence
    if sequences.keys() != coding:
        raise ValueError('Protein IDs do not match coding transcripts')
    if len(genes) != prediction['genes'] or len(coding) != prediction['coding_transcripts']:
        raise ValueError('Prediction counts differ from the accepted report')
    groups = defaultdict(list)
    for tid in sorted(coding):
        groups[parents[tid]].append(tid)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    selected = []
    with (output / f'{aid}.fa').open('w', encoding='utf-8') as fasta, \
         (output / 'representatives.tsv').open('w', encoding='utf-8') as table:
        table.write('gene_id\ttranscript_id\tprotein_length\tisoforms\n')
        for gid, tids in sorted(groups.items()):
            # Terminal stop is not an amino acid for length ranking. Retain the
            # original sequence, including any stops/partials, for tool assessment.
            tid = min(tids, key=lambda t: (-len(sequences[t].rstrip('*')), t))
            selected.append(tid)
            sequence = sequences[tid]
            fasta.write(f'>{tid}\n')
            for start in range(0, len(sequence), 60):
                fasta.write(sequence[start:start + 60] + '\n')
            table.write(f'{gid}\t{tid}\t{len(sequence.rstrip("*"))}\t{len(tids)}\n')
    data = dict(assembly_id=aid, status='ready', model_set='unfiltered',
                genes=len(genes), coding_genes=len(groups), coding_transcripts=len(coding),
                representative_proteins=len(selected),
                selection='Longest protein per coding gene; ties by transcript ID; terminal stop excluded from ranking',
                gff3_sha256=prediction['gff3_sha256'], proteins_sha256=prediction['proteins_sha256'],
                representative_sha256=sha(output / f'{aid}.fa'), settings=settings)
    (output / 'selection.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    return data


def report(selections, agat, busco, omark, output):
    indexed = {}
    for label, paths in [('agat', agat), ('busco', busco), ('omark', omark)]:
        rows = {}
        for item in paths:
            path = Path(item['path'])
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f'Empty/missing {label} result: {path}')
            aid = item['aid']
            if aid in rows:
                raise ValueError(f'Duplicate {label} result for {aid}')
            rows[aid] = path
        indexed[label] = rows
    result = []
    for path in selections:
        row = json.loads(Path(path).read_text(encoding='utf-8'))
        aid = row['assembly_id']
        if any(r['assembly_id'] == aid for r in result):
            raise ValueError('Duplicate baseline assembly')
        row['analyses'] = {}
        for label in ('agat', 'busco', 'omark'):
            enabled = label == 'agat' or row['settings'].get(f'run_{label}', False)
            if enabled and aid not in indexed[label]:
                raise ValueError(f'Missing requested {label} result for {aid}')
            row['analyses'][label] = 'ready' if enabled else 'not_requested'
        result.append(row)
    if not result:
        raise ValueError('No baseline assemblies')
    expected = {r['assembly_id'] for r in result}
    for label, rows in indexed.items():
        enabled = {r['assembly_id'] for r in result if r['analyses'][label] == 'ready'}
        if set(rows) != enabled or not set(rows) <= expected:
            raise ValueError(f'Unexpected {label} results')
    data = dict(schema_version=1, stage='annotation_baseline', status='ready',
                scope='Unfiltered coding models: AGAT all isoforms; BUSCO/OMArk one longest protein per coding gene; no model filtering',
                assemblies=sorted(result, key=lambda r: r['assembly_id']))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'baseline_status.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    lines = ['# Unfiltered annotation baseline', '', 'Status: **ready**', '', data['scope'], '',
             '| Assembly | Genes | Coding transcripts | Representative proteins | BUSCO | OMArk |',
             '|---|---:|---:|---:|---|---|']
    for row in data['assemblies']:
        aid = row['assembly_id']
        lines += [f"| {aid} | {row['genes']} | {row['coding_transcripts']} | {row['representative_proteins']} | {row['analyses']['busco']} | {row['analyses']['omark']} |", '',
                  f'[{aid}: AGAT statistics]({aid}/agat/{aid}_annotation_stats)',
                  f'[{aid}: representative selection]({aid}/selection/representatives.tsv)']
        if row['analyses']['busco'] == 'ready':
            lines.append(f'[{aid}: BUSCO summary]({aid}/busco/{aid}_busco.txt)')
        if row['analyses']['omark'] == 'ready':
            lines.append(f'[{aid}: OMArk summary]({aid}/omark/{aid}_detailed.txt)')
        lines.append('')
    lines += ['Representative selection does not remove or alter the original BRAKER models.',
              'Counts and isoform-sensitive completeness results must be interpreted with the stated protein set.',
              'BUSCO/OMArk metrics assess biological quality; ready means requested tools completed, not a quality cutoff.',
              'Protein-based AED and model filters are subsequent work; missing protein evidence is not evidence of an unsupported model.', '',
              '[Machine-readable status](baseline_status.json)']
    for label in ('busco', 'omark'):
        for aid, path in sorted(indexed[label].items()):
            lines += ['', f'## {aid}: {label.upper()} summary', '',
                      '```text', path.read_text(encoding='utf-8').strip(), '```']
    (output / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return data


def decode(value):
    return json.loads(base64.b64decode(value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'report'))
    parser.add_argument('--request', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    request = decode(args.request)
    if args.action == 'prepare':
        prepare(request['assembly_id'], request['settings'],
                json.loads(Path(request['validation']).read_text(encoding='utf-8')),
                request['gff'], request['proteins'], args.outdir, genome=request['genome'])
    else:
        report(request['selections'], request['agat'], request['busco'], request['omark'], args.outdir)


if __name__ == '__main__':
    main()
