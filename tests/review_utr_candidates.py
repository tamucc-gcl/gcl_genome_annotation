#!/usr/bin/env python3
"""Review completed UTR proposals without rerunning ranking or refinement."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from evidence_aed import merged, overlap_length, rna_hits
from inspect_utr_labels import annotations
from utr_refinement import subtract
from prepare_protein_evidence import checksum


def classify(old, new):
    """Separate serialization differences from changes to covered exon bases."""
    if sorted(old) == sorted(new):
        return 'identical_exon_records'
    if set(old) == set(new):
        return 'duplicate_exon_records_only'
    if merged(old) == merged(new):
        return 'same_exonic_bases_different_segmentation'
    added, removed = subtract(new, old), subtract(old, new)
    if added and removed:
        return 'added_and_removed_exonic_bases'
    return 'added_exonic_bases_only' if added else 'removed_exonic_bases_only'


def review(root, aid, output):
    root, output = Path(root), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    base = root / 'structural_annotation'
    folder = base / 'utr' / aid
    state = json.loads((folder / 'utr_status.json').read_text())
    identity = json.loads((root / 'filtered_assemblies' / aid / f'{aid}.contigs.json').read_text())
    lengths = {c['id']: c['length'] for c in identity['contigs']}
    if (state['status'] != 'ready' or state['reference_identity'] != identity['reference_identity'] or
        checksum(base / 'braker' / aid / f'{aid}.gff3') != state['original_model_sha256'] or
        checksum(folder / f'{aid}.utr.gff3') != state['gff3_sha256']):
        raise ValueError('Review inputs differ from the accepted UTR run')
    print('Reading original, native and accepted annotations...', flush=True)
    old = annotations(base / 'braker' / aid / f'{aid}.gff3', lengths)
    native = annotations(folder / 'native/native.gff3', lengths)
    accepted = annotations(folder / f'{aid}.utr.gff3', lengths)
    with (folder / 'utr_changes.tsv').open() as handle:
        changes = {r['transcript_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    rejected = state['rejected_utr_candidates']
    rows = []
    for tid, reason in rejected.items():
        a, b = old[tid], native[tid]
        rows.append(dict(transcript_id=tid, reason=reason,
            category=classify(a['exons'], b['exons']),
            added_bp=sum(e-s for s, e in subtract(b['exons'], a['exons'])),
            removed_bp=sum(e-s for s, e in subtract(a['exons'], b['exons'])),
            native_evidence_id=b['attributes'].get('utr_refine_evidence', '.'),
            original_cds_bp=sum(e-s for s,e in merged(a['cds'])),
            original_exonic_bp=sum(e-s for s,e in merged(a['exons'])),
            original_exon_records=len(a['exons']), native_exon_records=len(b['exons']),
            original_restored=(a['attributes'] == accepted[tid]['attributes'] and
                               sorted(a['exons']) == sorted(accepted[tid]['exons']))))
    with (output / 'rejected_proposals.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else
                                ['transcript_id', 'reason', 'category'], delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    selected = {}
    for col in ('five_prime_utr_bp', 'three_prime_utr_bp'):
        eligible = [r for r in changes.values() if r['changed'].lower() == 'true' and int(r[col]) > 0]
        for r in sorted(eligible, key=lambda r: (-int(r[col]), r['transcript_id']))[:5]:
            selected.setdefault(r['transcript_id'], []).append('long_' + col)
    for category in sorted({r['category'] for r in rows}):
        for r in [r for r in rows if r['category'] == category][:3]:
            selected.setdefault(r['transcript_id'], []).append('rejected_' + category)
    # Always include the known strand disagreement, even if another category dominates.
    for tid, reason in rejected.items():
        if 'contig or strand' in reason:
            selected.setdefault(tid, []).append('strand_disagreement')
    rna_ids = {accepted[t]['attributes'].get('utr_refine_evidence') for t in selected}
    rna_path = base / 'braker' / aid / 'braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
    rna = {}
    if rna_path.is_file():
        print('Reading selected RNA evidence...', flush=True)
        for hit in rna_hits(rna_path, lengths):
            if hit['identifier'] in rna_ids:
                rna[hit['identifier']] = hit
    # Group coding blocks by gene to distinguish span overlaps from coding-base overlaps.
    genes = {}
    for tid, model in old.items():
        gid = model['attributes']['Parent']
        g = genes.setdefault(gid, dict(seq=model['seq'], strand=model['strand'], cds=[]))
        g['cds'].extend(model['cds'])
    details = []
    for tid, why in selected.items():
        a, b, c = old[tid], native[tid], accepted[tid]
        extension = subtract(c['exons'], a['exons'])
        evidence_id = c['attributes'].get('utr_refine_evidence')
        hit = rna.get(evidence_id)
        neighbors = []
        lo, hi = min(s for s, e in c['exons']), max(e for s, e in c['exons'])
        for gid, gene in genes.items():
            if gid == a['attributes']['Parent'] or gene['seq'] != c['seq'] or not gene['cds']:
                continue
            if max(e for s, e in gene['cds']) <= lo or min(s for s, e in gene['cds']) >= hi:
                continue
            neighbors.append(dict(gene_id=gid, strand=gene['strand'],
                cds_span=[min(s for s, e in gene['cds']), max(e for s, e in gene['cds'])],
                extension_overlap_cds_bp=overlap_length(extension, gene['cds'])))
        def describe(model):
            return dict(seq=model['seq'], strand=model['strand'],
                feature_contexts=sorted(model['feature_contexts']),
                attributes=model['attributes'], exons=model['exons'], cds=model['cds'])
        details.append(dict(transcript_id=tid, selected_for=why, changes=changes[tid],
            rejection_reason=rejected.get(tid), original=describe(a), native=describe(b),
            accepted=describe(c), accepted_extension=extension,
            selected_rna=dict(id=evidence_id, seq=hit['seq'], strand=hit['strand'], blocks=hit['blocks']) if hit else None,
            extension_bases_outside_selected_rna=(sum(e-s for s,e in subtract(extension,hit['blocks'])) if hit else None),
            other_gene_coding_spans_overlapped=neighbors))
    result = dict(stage='utr_candidate_review', coordinates='0-based half-open',
        scope='Completed-file geometry and provenance review; BAM coverage not assessed; no annotation changes',
        transcripts=len(old), rejected=len(rows), rejected_categories=dict(Counter(r['category'] for r in rows)),
        all_rejected_originals_restored=all(r['original_restored'] for r in rows),
        retained_rna_found=rna_path.is_file(), selected_cases=details)
    (output / 'review.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'selected_cases'}, indent=2))
    print(f'Review files: {output}/review.json and rejected_proposals.tsv', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-root', required=True)
    parser.add_argument('--assembly-id', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    review(args.inputs_root, args.assembly_id, args.outdir)
