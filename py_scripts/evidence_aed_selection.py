#!/usr/bin/env python3
"""AED-style intact-chain selection, with the original broad union as a control.

This assesses the closest available RNA/protein structures; it is not independent
validation or a reimplementation of MAKER's complete evidence selection.
Normalization and expensive alignments remain unchanged.
"""
import argparse
from contextlib import closing
import csv
import json
from pathlib import Path
import sqlite3
import statistics

from evidence_aed import MAKER_REVISION, models, merged, nucleotide_aed, overlap_length
from prepare_protein_evidence import checksum


def select_chains(model_blocks, candidates):
    """One intact best-matching chain per type; do not concatenate isoforms.

    Deduplicate identical block structures before ranking. AED uses full blocks,
    including overhangs; no clipping to the annotation. Ties use source ID then
    alignment ID for deterministic choice, not database iteration order.
    """
    unique = {}
    for hit in candidates:
        blocks = tuple(merged(hit['blocks']))
        if not overlap_length(model_blocks, blocks):
            continue
        hit = dict(hit, blocks=blocks)
        hit['metrics'] = nucleotide_aed(model_blocks, blocks)
        previous = unique.get(blocks)
        if previous is None or (hit['target'], hit['identifier']) < (previous['target'], previous['identifier']):
            unique[blocks] = hit
    ordered = sorted(unique.values(), key=lambda h: (h['metrics']['aed'], h['target'], h['identifier']))
    return (ordered[0] if ordered else None), len(unique)


def score(gff, identity, normalized, output, validation=None, selection='best_chain'):
    if selection not in ('best_chain', 'all_union'):
        raise ValueError('Unknown AED evidence selection policy')
    normalized, output = Path(normalized), Path(output)
    metadata = json.loads((normalized / 'evidence.json').read_text(encoding='utf-8'))
    if metadata['reference_identity'] != identity['reference_identity'] or checksum(normalized / 'evidence.sqlite') != metadata['database_sha256']:
        raise ValueError('Evidence reference/checksum mismatch')
    if validation is not None and (validation.get('status') != 'ready' or
            validation.get('gff3_sha256') != checksum(gff) or
            validation.get('reference_identity') != identity['reference_identity']):
        raise ValueError('Models/reference differ from accepted BRAKER validation')
    transcripts = models(gff, {c['id']: c['length'] for c in identity['contigs']})
    available = {k: v['available'] for k, v in metadata['evidence'].items()}
    rows = []
    with closing(sqlite3.connect(f'{(normalized / "evidence.sqlite").resolve().as_uri()}?mode=ro', uri=True)) as connection:
        for tid, model in sorted(transcripts.items()):
            seq = metadata['coordinates'][model['seq']]
            strand = 1 if model['strand'] == '+' else -1
            candidates = {'rna': [], 'protein': []}
            query = connection.execute('SELECT kind,identifier,target,blocks FROM spans JOIN hits USING(id) '
                'WHERE seq_min<=? AND seq_max>=? AND strand_min<=? AND strand_max>=? AND start<? AND end>?',
                (seq, seq, strand, strand, model['end'], model['start']))
            for kind, identifier, target, text in query:
                blocks = json.loads(text)
                # Remove RTree floating-point boundary false positives exactly.
                if max(e for s, e in blocks) <= model['start'] or min(s for s, e in blocks) >= model['end']:
                    continue
                candidates[kind].append(dict(identifier=identifier, target=target, blocks=blocks))
            best, distinct, broad, selected = {}, {}, {}, {}
            for kind in ('rna', 'protein'):
                best[kind], distinct[kind] = select_chains(model['exons'], candidates[kind])
                broad[kind] = merged([b for h in candidates[kind] for b in h['blocks']])
                selected[kind] = (list(best[kind]['blocks']) if best[kind] else []) if selection == 'best_chain' else broad[kind]
            result = nucleotide_aed(model['exons'], selected['rna'] + selected['protein'])
            original = nucleotide_aed(model['exons'], broad['rna'] + broad['protein'])
            row = dict(gene_id=model['gene_id'], transcript_id=tid,
                aed=result['aed'] if any(available.values()) else None,
                combined_aed=result['aed'] if any(available.values()) else None,
                status='scored' if any(available.values()) else 'unavailable',
                **{k: result[k] for k in ('sensitivity','specificity','model_bp','evidence_bp','overlap_bp')},
                rna_hits=len(candidates['rna']), protein_hits=len(candidates['protein']))
            for kind in ('rna', 'protein'):
                row[f'{kind}_aed'] = nucleotide_aed(model['exons'], selected[kind])['aed'] if available[kind] else None
            row['legacy_protein_cds_aed'] = (min([1.0] + [nucleotide_aed(model['cds'], h['blocks'])['aed']
                for h in candidates['protein']]) if available['protein'] and model['cds'] else None)
            row['all_union_combined_aed'] = original['aed'] if any(available.values()) else None
            row['all_union_evidence_bp'] = original['evidence_bp']
            for kind in ('rna', 'protein'):
                row[f'all_union_{kind}_aed'] = nucleotide_aed(model['exons'], broad[kind])['aed'] if available[kind] else None
                hit = best[kind]
                row[f'{kind}_distinct_overlapping_chains'] = distinct[kind]
                row[f'{kind}_best_alignment_id'] = hit['identifier'] if hit else None
                row[f'{kind}_best_source_id'] = hit['target'] if hit else None
                row[f'{kind}_selected_hits'] = int(hit is not None) if selection == 'best_chain' else len(candidates[kind])
            row['selection_policy'] = selection
            row['reference_contig'] = model['seq']
            row['strand'] = model['strand']
            row['model_exon_blocks'] = json.dumps(merged(model['exons']), separators=(',', ':'))
            row['model_cds_blocks'] = json.dumps(merged(model['cds']), separators=(',', ':'))
            row['model_cds_bp'] = sum(e-s for s,e in merged(model['cds']))
            row['model_exon_outside_cds_bp'] = row['model_bp'] - overlap_length(model['exons'], model['cds'])
            for kind in ('rna', 'protein'):
                blocks = selected[kind]
                metrics = nucleotide_aed(model['exons'], blocks)
                row[f'{kind}_selected_blocks'] = json.dumps(blocks, separators=(',', ':'))
                row[f'{kind}_model_coverage'] = metrics['specificity'] if available[kind] else None
                row[f'{kind}_evidence_outside_model_bp'] = metrics['evidence_bp']-metrics['overlap_bp'] if available[kind] else None
            rows.append(row)
    output.mkdir(parents=True, exist_ok=True)
    primary = ['gene_id', 'transcript_id', 'rna_aed', 'protein_aed', 'combined_aed', 'legacy_protein_cds_aed',
        'status', 'sensitivity', 'specificity', 'model_bp', 'evidence_bp', 'overlap_bp',
        'rna_hits', 'protein_hits', 'rna_distinct_overlapping_chains',
        'rna_best_alignment_id', 'rna_best_source_id', 'rna_selected_hits',
        'protein_distinct_overlapping_chains', 'protein_best_alignment_id',
        'protein_best_source_id', 'protein_selected_hits', 'selection_policy',
        'model_cds_bp', 'model_exon_outside_cds_bp', 'rna_model_coverage',
        'rna_evidence_outside_model_bp', 'protein_model_coverage', 'protein_evidence_outside_model_bp']
    diagnostics = ['gene_id', 'transcript_id'] + [k for k in rows[0] if k not in primary and k != 'aed']
    for filename, fields in [('aed_scores.tsv', primary), ('aed_diagnostics.tsv', diagnostics)]:
        with (output / filename).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            for row in rows:
                writer.writerow({k: '.' if row[k] is None else f'{row[k]:.8f}' if isinstance(row[k], float) else row[k] for k in fields})
    genes = {}
    for row in rows:
        genes.setdefault(row['gene_id'], []).append(row['aed'])
    summary = dict(schema_version=3, status='ready', transcripts=len(rows), genes=len(genes),
        maker_revision=MAKER_REVISION, selection_policy=selection,
        definition='AED nucleotide arithmetic on intact selected evidence chains; not full MAKER AED',
        selection='Best individual full RNA and protein chain by transcript-exon AED; union these two for combined score' if selection == 'best_chain' else 'Original full union of same-contig/strand span-overlapping chains',
        limitations=['Selection is annotation-conditioned and optimistic; not a probability or independent validation',
                    'Single-chain selection can miss complementary fragments within the same evidence type',
                    'Best RNA and protein chains can conflict; combined score need not improve',
                    'Exon-level overlap does not establish splice-site, reading-frame or sequence correctness'],
        evidence=metadata['evidence'], model_sha256=checksum(gff), reference_identity=identity['reference_identity'],
        gene_aed={g:min(v for v in values if v is not None) if any(v is not None for v in values) else None for g,values in genes.items()},
        comparison=dict(combined_worse_than_both=sum(r['combined_aed'] > max(r['rna_aed'],r['protein_aed']) + 1e-7
            for r in rows if all(r[k] is not None for k in ('combined_aed','rna_aed','protein_aed'))),
            full_model_coverage=sum(r['model_bp']==r['overlap_bp'] for r in rows)))
    summary['score_distributions'] = {}
    for key in ('rna_aed','protein_aed','combined_aed'):
        values = [r[key] for r in rows if r[key] is not None]
        summary['score_distributions'][key] = dict(available=len(values),
            median=statistics.median(values) if values else None,
            mean=statistics.mean(values) if values else None,
            above_0_7=sum(v>.7 for v in values), no_overlap=sum(v==1 for v in values))
    (output / 'aed_status.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    (output / 'README.md').write_text(
        '# Evidence agreement assessment\n\nStatus: **ready**\n\n'
        f'Selection policy: `{selection}`. Models are unchanged; no filters are applied.\n\n'
        '[Transcript scores](aed_scores.tsv) | [Definition and provenance](aed_status.json) | [Evaluation diagnostics](aed_diagnostics.tsv)\n\n'
        '`rna_aed`, `protein_aed`, and `combined_aed` are AED-style agreement scores (lower is better). '
        'The score table retains supporting coverage, overlap, hit counts and selected evidence IDs. '
        '`legacy_protein_cds_aed` is temporarily retained in the main table for validation; the duplicate `aed` alias is omitted. '
        'Missing sources are `.`, supplied but unsupported sources score 1.\n\n'
        'With best_chain, select the closest intact RNA and protein alignment separately, then union those two for combined AED. '
        'Full blocks remain, including overhangs and conflicting exons; blocks are never clipped to the annotation. '
        'Identical chains are deduplicated for selection. Intronic span-only matches cannot supply exon support.\n\n'
        'The separate diagnostics table retains historical comparisons and selected evidence/model coordinate blocks '
        'for evaluation. Coordinates are 0-based half-open. '
        'Historical scores are not additional production scores.\n\n'
        'This asks whether a matching evidence chain exists and can be optimistic. It can miss complementary '
        'fragments within one source type; combined evidence can still disagree. These scores are not full MAKER AED, '
        'a probability of correctness, or independent validation. Splice junction assessment follows separately.\n', encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('models','identity','normalized','validation','outdir'):
        parser.add_argument('--'+arg, required=True)
    args = parser.parse_args()
    score(args.models, json.loads(Path(args.identity).read_text()), args.normalized, args.outdir,
          json.loads(Path(args.validation).read_text()))


if __name__ == '__main__':
    main()
