#!/usr/bin/env python3
"""Validate inputs and report native InGenAnnot scores; no custom AED arithmetic."""
import argparse
import csv
import json
from pathlib import Path
import statistics

from evidence_aed import features, models
from prepare_protein_evidence import checksum
from ingenannot_inputs import check


def report(original, annotated, identity, inputs, current_scores, cds_only, penalty, output):
    identity = json.loads(Path(identity).read_text())
    accepted = json.loads(Path(inputs).read_text())
    if accepted.get('status') != 'ready' or identity['reference_identity'] != accepted['reference_identity'] or checksum(original) != accepted['model_sha256']:
        raise ValueError('Benchmark reference/models differ from validated inputs')
    lengths = {c['id']:c['length'] for c in identity['contigs']}
    before, after = models(original,lengths), models(annotated,lengths)
    if set(before) != set(after):
        raise ValueError('InGenAnnot changed transcript identities during assessment')
    normalized_bounds = []
    for tid,a in before.items():
        b = after[tid]
        if any(a[k]!=b[k] for k in ('gene_id','seq','strand')) or any(sorted(a[k])!=sorted(b[k]) for k in ('exons','cds')):
            raise ValueError(f'InGenAnnot changed model structure during assessment: {tid}')
        original_bounds = (a['start'],a['end'])
        native_bounds = (b['start'],b['end'])
        if original_bounds != native_bounds:
            # InGenAnnot 0.0.20 Utils fixes transcript bounds to the exon envelope.
            # BRAKER may include a separately recorded stop codon in its mRNA span.
            exon_bounds = (min(s for s,e in a['exons']),max(e for s,e in a['exons']))
            if native_bounds != exon_bounds:
                raise ValueError(f'InGenAnnot changed model structure during assessment: {tid}: unexpected transcript bounds')
            normalized_bounds.append(tid)
    def coding_phases(path):
        result = {tid:[] for tid in before}
        for row in features(path,lengths):
            if row['kind'] == 'CDS':
                for tid in row['attrs']['Parent'].split(','):
                    result[tid].append((row['start'],row['end'],row['phase']))
        return {tid:sorted(parts) for tid,parts in result.items()}
    if coding_phases(original) != coding_phases(annotated):
        raise ValueError('InGenAnnot changed model structure during assessment: coding phases differ')
    scores = {}
    for row in features(annotated,lengths):
        if row['kind'] not in ('mRNA','transcript'): continue
        attrs = row['attrs']
        tid = attrs['ID']
        entry = {}
        for kind,field,evidence in [('rna','aed_ev_tr','ev_tr'),('protein','aed_ev_pr','ev_pr')]:
            value = float(attrs[field])
            if not 0 <= value <= 1:
                raise ValueError('Invalid native AED value')
            entry[kind] = value if accepted['evidence'][kind]['available'] else None
            entry[kind+'_evidence_id'] = attrs.get(evidence)
        entry['rna_penalty'] = attrs.get('ev_tr_penalty')
        scores[tid] = entry
    with Path(current_scores).open() as handle:
        baseline = {r['transcript_id']:r for r in csv.DictReader(handle,delimiter='\t')}
    if set(baseline) != set(before) or set(scores) != set(before):
        raise ValueError('Benchmark score IDs differ from model IDs')
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    rows = []
    for tid in sorted(before):
        native = scores[tid]
        row = dict(gene_id=before[tid]['gene_id'],transcript_id=tid)
        if baseline[tid]['gene_id'] != row['gene_id']:
            raise ValueError('Benchmark gene ID differs from current AED')
        for kind in ('rna','protein'):
            old = baseline[tid][kind+'_aed']
            old = None if old=='.' else float(old)
            if (old is not None) != accepted['evidence'][kind]['available']:
                raise ValueError('Baseline source availability mismatch')
            row['current_'+kind+'_aed'] = old
            row['ingenannot_'+kind+'_aed'] = native[kind]
            row[kind+'_delta'] = native[kind]-old if old is not None else None
            row[kind+'_evidence_id'] = native[kind+'_evidence_id']
        row['rna_penalty'] = native['rna_penalty']
        rows.append(row)
    with (output/'benchmark_scores.tsv').open('w',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t'); writer.writeheader()
        for row in rows:
            writer.writerow({k:'.' if v is None else f'{v:.8f}' if isinstance(v,float) else v for k,v in row.items()})
    summary = dict(status='ready',stage='ingenannot_benchmark',transcripts=len(rows),genes=len({r['gene_id'] for r in rows}),
        models_unchanged=True,transcript_bounds_normalized_to_exons=normalized_bounds,
        model_validation_scope='Transcript IDs, gene parents, reference, strand, exact exon/CDS blocks and CDS phases preserved; native transcript bounds may normalize to the exon envelope',
        profile=dict(version='0.0.20',rna_cds_only=cds_only,rna_overflow_penalty=penalty,protein_cds_only=True),
        model_sha256=checksum(original),reference_identity=identity['reference_identity'],evidence=accepted['evidence'],metrics={})
    for key in ('current_rna_aed','ingenannot_rna_aed','current_protein_aed','ingenannot_protein_aed'):
        values = [r[key] for r in rows if r[key] is not None]
        summary['metrics'][key] = dict(available=len(values),median=statistics.median(values) if values else None,
            mean=statistics.mean(values) if values else None,score_equals_one=sum(v==1 for v in values),above_0_7=sum(v>.7 for v in values))
        if key.startswith('ingenannot_'):
            kind=key.split('_')[1]
            summary['metrics'][key]['no_selected_evidence']=sum(
                r[key] is not None and r[kind+'_evidence_id'] in (None,'None','.') for r in rows)
            summary['metrics'][key]['score_one_with_previous_overlap']=sum(
                r[key]==1 and r['current_'+kind+'_aed']<1 for r in rows if r[key] is not None)
    (output/'benchmark_status.json').write_text(json.dumps(summary,indent=2)+'\n')
    (output/'README.md').write_text('# InGenAnnot AED benchmark\n\nStatus: **ready**\n\n'
        '[Matched RNA/protein scores](benchmark_scores.tsv) | [Settings and provenance](benchmark_status.json)\n\n'
        'Native InGenAnnot assesses unchanged BRAKER models using the same RNA/protein alignments as the current method. '
        'No combined AED is produced. Differences include model scope, evidence handling and the configured RNA penalty; '
        'lower scores alone do not establish improved accuracy. Native GFF scores are rounded to four decimal places.\n\n'
        f'InGenAnnot normalized {len(normalized_bounds)} transcript bounds to their unchanged exon envelopes. '
        'Transcript IDs, gene parents, reference, strand, exon/CDS blocks and coding phases are preserved. '
        'The original BRAKER GFF remains the annotation; the native GFF is an assessment export.\n\n'
        'The paper-based profile uses CDS-only RNA assessment and a 0.25 overflow/splice-disagreement penalty; '
        'these choices are benchmark settings, not calibrated fish filtering thresholds. No models are filtered or refined.\n\n'
        'A native score of 1 does not establish zero coordinate overlap. Penalized candidates can fail to beat '
        'the initial score of 1 and leave no evidence selected. Summary counts distinguish score 1, missing '
        'selected evidence and score-1 models that had overlap in the unpenalized baseline.\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action',required=True)
    p = sub.add_parser('check')
    for arg in ('gff','validation','current-status','out'): p.add_argument('--'+arg,required=True)
    p.add_argument('--rna'); p.add_argument('--protein')
    p = sub.add_parser('report')
    for arg in ('original','annotated','identity','inputs','current-scores','outdir'): p.add_argument('--'+arg,required=True)
    p.add_argument('--cds-only',choices=('true','false'),required=True)
    p.add_argument('--penalty',type=float,required=True)
    a = parser.parse_args()
    if a.action=='check': check(a.gff,a.validation,a.current_status,a.rna,a.protein,a.out)
    else: report(a.original,a.annotated,a.identity,a.inputs,a.current_scores,a.cds_only=='true',a.penalty,a.outdir)


if __name__=='__main__': main()
