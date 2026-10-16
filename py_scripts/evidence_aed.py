#!/usr/bin/env python3
"""Nucleotide AED following MAKER shadow_AED::get_AED/get_SN_SP.

Evidence selection is explicit: same reference and strand, overlapping the
supplied model span. No MAKER prediction, polishing, ORF merging or evidence
cluster heuristics are applied. Given identical selected blocks the core score
matches the pinned routine; this is not a claim of whole-pipeline equivalence.
"""
import argparse
from contextlib import closing
import csv
import json
from pathlib import Path
import re
import sqlite3
from urllib.parse import unquote

from prepare_protein_evidence import checksum

MAKER_REVISION = '4bf3085ae8df4680468646a63c635094a0a3dc2a'


def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError('Invalid half-open interval')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def overlap_length(left, right):
    left, right = merged(left), merged(right)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def nucleotide_aed(model_blocks, evidence_blocks):
    """MAKER's union/intersection arithmetic without chromosome-sized vectors."""
    model, evidence = merged(model_blocks), merged(evidence_blocks)
    if not model:
        raise ValueError('Cannot score an empty model')
    predicted = sum(e - s for s, e in model)
    supported = sum(e - s for s, e in evidence)
    common = overlap_length(model, evidence)
    sensitivity = common / supported if supported else 0.0
    specificity = common / predicted
    return dict(aed=1 - (sensitivity + specificity) / 2,
                sensitivity=sensitivity, specificity=specificity,
                model_bp=predicted, evidence_bp=supported, overlap_bp=common)


def features(path, lengths, gtf=False):
    with Path(path).open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            if line.startswith('>'):
                raise ValueError('Embedded FASTA is not supported in evidence GFF')
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected 9 fields')
            seq, _, kind, begin, finish, _, strand, phase, text = fields
            start, end = int(begin) - 1, int(finish)
            if seq not in lengths or not 0 <= start < end <= lengths[seq]:
                raise ValueError(f'{path}:{number}: feature outside accepted reference')
            if strand not in (('+', '-', '.') if gtf else ('+', '-')):
                raise ValueError(f'{path}:{number}: unresolved feature strand')
            attrs = (dict(re.findall(r'(\w+)\s+"([^"\n]*)"', text)) if gtf else
                     {k: unquote(v) for k, v in (x.split('=', 1) for x in text.split(';') if '=' in x)})
            yield dict(seq=seq, kind=kind, start=start, end=end, strand=strand, phase=phase, attrs=attrs)


def models(path, lengths):
    genes, transcripts, parts = {}, {}, []
    for row in features(path, lengths):
        kind, attrs = row['kind'], row['attrs']
        identifier = attrs.get('ID')
        parents = attrs.get('Parent', '').split(',')
        if kind == 'gene':
            if not identifier or identifier in genes:
                raise ValueError('Duplicate/missing model gene ID')
            genes[identifier] = row
        elif kind in ('mRNA', 'transcript'):
            if not identifier or identifier in transcripts or len(parents) != 1 or not parents[0]:
                raise ValueError('Model transcript requires unique ID and explicit gene parent')
            transcripts[identifier] = row | dict(gene_id=parents[0], exons=[], cds=[])
        elif kind in ('exon', 'CDS'):
            if not parents[0]:
                raise ValueError('Model block missing Parent')
            parts.append((row, parents))
    for row, parents in parts:
        for parent in parents:
            if parent not in transcripts:
                raise ValueError('Unresolved model block Parent')
            model = transcripts[parent]
            if (row['seq'], row['strand']) != (model['seq'], model['strand']) or not model['start'] <= row['start'] < row['end'] <= model['end']:
                raise ValueError(f"Model block differs from transcript coordinates/strand: {parent}; "
                    f"{row['kind']} {row['seq']}:{row['start']+1}-{row['end']} ({row['strand']}); "
                    f"transcript {model['seq']}:{model['start']+1}-{model['end']} ({model['strand']})")
            model['exons' if row['kind'] == 'exon' else 'cds'].append((row['start'], row['end']))
    for tid, row in transcripts.items():
        if row['gene_id'] not in genes or not row['exons']:
            raise ValueError(f'{tid}: unresolved gene or missing transcript exons')
        gene = genes[row['gene_id']]
        if (row['seq'], row['strand']) != (gene['seq'], gene['strand']) or not gene['start'] <= row['start'] < row['end'] <= gene['end']:
            raise ValueError('Transcript differs from gene bounds/strand')
    if not transcripts:
        raise ValueError('No models to score')
    return transcripts


def rna_hits(path, lengths):
    transcripts, exons = {}, []
    for row in features(path, lengths, gtf=True):
        tid = row['attrs'].get('transcript_id')
        if row['kind'] == 'transcript':
            if not tid or tid in transcripts:
                raise ValueError('Duplicate/missing RNA transcript ID')
            transcripts[tid] = row | dict(identifier=tid, blocks=[])
        elif row['kind'] == 'exon':
            exons.append((tid, row))
    for tid, exon in exons:
        if tid not in transcripts:
            raise ValueError('RNA exon has unresolved transcript')
        hit = transcripts[tid]
        if (exon['seq'], exon['strand']) != (hit['seq'], hit['strand']) or not hit['start'] <= exon['start'] < exon['end'] <= hit['end']:
            raise ValueError('RNA exon outside transcript bounds/strand')
        hit['blocks'].append((exon['start'], exon['end']))
    for hit in transcripts.values():
        if not hit['blocks']:
            raise ValueError('RNA transcript has no exons')
        yield hit


def protein_hits(path, lengths):
    """Stream miniprot's mRNA/CDS groups without retaining the whole alignment."""
    current = None
    for row in features(path, lengths):
        if row['kind'] == 'mRNA':
            if current:
                if not current['blocks']:
                    raise ValueError('Protein alignment has no CDS blocks')
                yield current
            identifier = row['attrs'].get('ID')
            target = row['attrs'].get('Target', '').split()
            if not identifier or not target:
                raise ValueError('Protein alignment missing ID/Target')
            current = row | dict(identifier=identifier, target=target[0], blocks=[])
        elif row['kind'] == 'CDS':
            if not current or row['attrs'].get('Parent') != current['identifier']:
                raise ValueError('Protein CDS is not grouped with its mRNA')
            if (row['seq'], row['strand']) != (current['seq'], current['strand']) or not current['start'] <= row['start'] < row['end'] <= current['end']:
                raise ValueError('Protein CDS outside alignment bounds/strand')
            current['blocks'].append((row['start'], row['end']))
    if current:
        if not current['blocks']:
            raise ValueError('Protein alignment has no CDS blocks')
        yield current


def normalize(identity, rna, protein, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lengths = {c['id']: c['length'] for c in identity['contigs']}
    coordinates = {name: i for i, name in enumerate(sorted(lengths))}
    database = output / 'evidence.sqlite'
    if database.exists():
        raise ValueError('Refusing to overwrite normalized evidence')
    counts = {}
    with closing(sqlite3.connect(database)) as connection:
        connection.executescript('CREATE TABLE hits(id INTEGER PRIMARY KEY, kind TEXT, identifier TEXT, target TEXT, blocks TEXT);'
                                'CREATE VIRTUAL TABLE spans USING rtree(id, seq_min, seq_max, strand_min, strand_max, start, end);')
        for kind, path, reader in [('rna', rna, rna_hits), ('protein', protein, protein_hits)]:
            count = ambiguous = 0
            if path is not None:
                for hit in reader(path, lengths):
                    if hit['strand'] == '.':
                        ambiguous += 1
                        continue
                    blocks = merged(hit['blocks'])
                    cursor = connection.execute('INSERT INTO hits(kind,identifier,target,blocks) VALUES(?,?,?,?)',
                                                (kind, hit['identifier'], hit.get('target', hit['identifier']), json.dumps(blocks)))
                    strand = 1 if hit['strand'] == '+' else -1
                    seq = coordinates[hit['seq']]
                    connection.execute('INSERT INTO spans VALUES(?,?,?,?,?,?,?)',
                                       (cursor.lastrowid, seq, seq, strand, strand, hit['start'], hit['end']))
                    count += 1
            counts[kind] = dict(available=path is not None, hits=count, ambiguous_strand_excluded=ambiguous,
                                source_sha256=checksum(path) if path is not None else None)
        connection.commit()
    data = dict(schema_version=1, reference_identity=identity['reference_identity'], coordinates=coordinates,
                evidence=counts, database_sha256=checksum(database))
    (output / 'evidence.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    return data


def score(gff, identity, normalized, output, validation=None):
    normalized, output = Path(normalized), Path(output)
    metadata = json.loads((normalized / 'evidence.json').read_text(encoding='utf-8'))
    if metadata['reference_identity'] != identity['reference_identity'] or checksum(normalized / 'evidence.sqlite') != metadata['database_sha256']:
        raise ValueError('Evidence reference/checksum mismatch')
    lengths = {c['id']: c['length'] for c in identity['contigs']}
    transcripts = models(gff, lengths)
    if validation is not None:
        if validation.get('status') != 'ready' or validation.get('gff3_sha256') != checksum(gff) or validation.get('reference_identity') != identity['reference_identity']:
            raise ValueError('Models/reference differ from accepted BRAKER validation')
    output.mkdir(parents=True, exist_ok=True)
    available = {kind: row['available'] for kind, row in metadata['evidence'].items()}
    rows = []
    with closing(sqlite3.connect(f'{(normalized / "evidence.sqlite").resolve().as_uri()}?mode=ro', uri=True)) as connection:
        for tid, model in sorted(transcripts.items()):
            seq = metadata['coordinates'][model['seq']]
            strand = 1 if model['strand'] == '+' else -1
            candidates = connection.execute('SELECT kind,identifier,blocks FROM spans JOIN hits USING(id) '
                                            'WHERE seq_min<=? AND seq_max>=? AND strand_min<=? AND strand_max>=? AND start<? AND end>?',
                                            (seq, seq, strand, strand, model['end'], model['start']))
            selected = {'rna': [], 'protein': []}
            hits = {'rna': 0, 'protein': 0}
            for kind, identifier, text in candidates:
                blocks = json.loads(text)
                # RTree stores rounded float bounds; enforce exact overlap.
                if max(e for s, e in blocks) <= model['start'] or min(s for s, e in blocks) >= model['end']:
                    continue
                selected[kind].extend(blocks)
                hits[kind] += 1
            result = nucleotide_aed(model['exons'], selected['rna'] + selected['protein'])
            row = dict(gene_id=model['gene_id'], transcript_id=tid,
                       aed=result['aed'] if any(available.values()) else None,
                       combined_aed=result['aed'] if any(available.values()) else None,
                       status='scored' if any(available.values()) else 'unavailable',
                       **{k: result[k] for k in ('sensitivity', 'specificity', 'model_bp', 'evidence_bp', 'overlap_bp')},
                       rna_hits=hits['rna'], protein_hits=hits['protein'])
            for kind in ('rna', 'protein'):
                row[f'{kind}_aed'] = nucleotide_aed(model['exons'], selected[kind])['aed'] if available[kind] else None
            row['legacy_protein_cds_aed'] = None
            # Legacy comparison uses the best individual protein alignment, not
            # the canonical union, and CDS rather than transcript exon blocks.
            if available['protein'] and model['cds']:
                proteins = connection.execute('SELECT blocks FROM spans JOIN hits USING(id) '
                                              'WHERE kind="protein" AND seq_min<=? AND seq_max>=? AND strand_min<=? AND strand_max>=? AND start<? AND end>?',
                                              (seq, seq, strand, strand, model['end'], model['start']))
                row['legacy_protein_cds_aed'] = min([1.0] + [nucleotide_aed(model['cds'], json.loads(p[0]))['aed'] for p in proteins])
            rows.append(row)
    with (output / 'aed_scores.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: '.' if v is None else f'{v:.8f}' if isinstance(v, float) else v for k, v in row.items()})
    genes = {}
    for row in rows:
        genes.setdefault(row['gene_id'], []).append(row['aed'])
    summary = dict(schema_version=1, status='ready', transcripts=len(rows), genes=len(genes),
                   maker_revision=MAKER_REVISION, definition='shadow_AED::get_AED/get_SN_SP: union of selected blocks',
                   selection='Same contig/strand, evidence span overlaps each unchanged transcript; full blocks retained including outside model span',
                   compatibility_scope='MAKER nucleotide core on identical selected blocks; full MAKER evidence clustering, polishing, eAED and QI not replicated',
                   evidence=metadata['evidence'], model_sha256=checksum(gff), reference_identity=identity['reference_identity'],
                   gene_aed={gid: min(v for v in values if v is not None) if any(v is not None for v in values) else None for gid, values in genes.items()})
    (output / 'aed_status.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    (output / 'README.md').write_text(
        '# Evidence AED assessment\n\nStatus: **ready**\n\n'
        f"Scored {len(rows)} transcripts from {len(genes)} genes. Gene summaries use the lowest available transcript AED.\n\n"
        '[Transcript scores](aed_scores.tsv) | [Scoring definition and provenance](aed_status.json)\n\n'
        'The three evidence-agreement scores are `rna_aed`, `protein_aed` and `combined_aed`. '
        'Lower values indicate greater agreement (0 = complete agreement; 1 = no overlap). '
        '`aed` is retained as an alias of `combined_aed`. These are AED-style agreement scores, not probabilities of correctness.\n\n'
        'AED = 1 − (sensitivity + specificity) / 2. RNA exon blocks and protein alignment blocks are unioned before scoring; repeated hits do not add weight. '
        'Models are unchanged. Scores use transcript exons, including annotated UTRs, and full selected evidence blocks.\n\n'
        'Evidence selection requires matching contig and strand and overlap with the transcript span. '
        'The nucleotide arithmetic follows MAKER; MAKER’s full evidence clustering, eAED and annotation polishing are not reproduced.\n\n'
        'A missing evidence source is unavailable (`.`), while supplied evidence with no matching support gives AED 1. '
        'RNA transcripts with ambiguous strand are excluded and counted in the status report. '
        'The legacy protein/CDS comparison is reported separately. No filters are applied. '
        'Agreement with training evidence is not independent biological validation.\n', encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('normalize', 'score'))
    parser.add_argument('--identity', required=True)
    parser.add_argument('--rna')
    parser.add_argument('--protein')
    parser.add_argument('--models')
    parser.add_argument('--normalized')
    parser.add_argument('--validation')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    identity = json.loads(Path(args.identity).read_text(encoding='utf-8'))
    if args.action == 'normalize':
        normalize(identity, args.rna, args.protein, args.outdir)
    else:
        score(args.models, identity, args.normalized, args.outdir,
              json.loads(Path(args.validation).read_text(encoding='utf-8')) if args.validation else None)


if __name__ == '__main__':
    main()
