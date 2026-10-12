#!/usr/bin/env python3
"""Prediction summary separate from the cached BRAKER training task."""
import argparse
import json
from pathlib import Path


def summarize(plan, predictions, runtimes, check_only, outdir, rna_qc_link='../rna_alignment/multiqc/multiqc_report.html'):
    results = {r['assembly_id']: r for r in predictions}
    checks = {r['assembly_id']: r for r in runtimes}
    rows = []
    for job in plan['jobs']:
        aid = job['assembly_id']
        row = {'assembly_id': aid, 'mode': job['mode'], 'readsets': job['readsets'],
               'protein_count': job['protein_count'], 'reference_identity': job['reference_identity']}
        if aid not in checks:
            row.update(status='blocked', reason='BRAKER software check missing')
        elif check_only:
            row.update(status='software_checked', reason='Training and prediction were not requested')
        elif aid not in results or results[aid]['status'] != 'ready':
            row.update(status='blocked', reason='Validated BRAKER predictions missing')
        else:
            row.update(results[aid])
        rows.append(row)
    status = 'incomplete' if any(r['status'] == 'blocked' for r in rows) else ('software_checked' if check_only else 'ready')
    result = {'schema_version': 1, 'stage': 'braker_prediction', 'status': status, 'assemblies': rows,
              'scope': 'BRAKER3 coding-gene prediction; downstream filtering, functional annotation and biological quality assessment remain pending'}
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    (outdir / 'prediction_status.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    lines = ['# BRAKER3 prediction', '', f'Status: **{status}**', '',
             '| Assembly | Status | Genes | Coding transcripts |', '|---|---|---:|---:|']
    for row in rows:
        aid = row['assembly_id']
        lines.append(f'| {aid} | {row["status"]} | {row.get("genes", "pending")} | {row.get("coding_transcripts", "pending")} |')
        lines.extend(['', f'[{aid}: runtime log](braker_runtime/{aid}/console.log)',
                      f'[{aid}: validated input handoff](braker_inputs/{aid}.inputs.json)'])
        if row['status'] == 'ready':
            lines.extend([f'[{aid}: gene models](braker/{aid}/{aid}.gff3)',
                          f'[{aid}: proteins](braker/{aid}/{aid}_proteins.fa)',
                          f'[{aid}: CDS](braker/{aid}/{aid}_cds.fa)',
                          f'[{aid}: execution and training log](braker/{aid}/console.log)',
                          f'[{aid}: command and evidence provenance](braker/{aid}/inputs.json)'])
        if row.get('reason'):
            lines.extend(['', row['reason']])
    lines.extend(['', '[Resolved prediction status](prediction_status.json)',
                  f'[RNA alignment QC]({rna_qc_link})', '',
                  'Readset BAMs remain separate. BRAKER3 handles transcript assembly internally.',
                  'Library strand metadata is retained; UTR prediction and strand-separated BAM flags are not enabled.',
                  'These are BRAKER outputs before the pipeline model filters. Structural and functional evaluation follow.', ''])
    (outdir / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('plan', 'predictions', 'runtimes', 'check-only', 'outdir'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--rna-qc-link', default='../rna_alignment/multiqc/multiqc_report.html')
    args = parser.parse_args()
    load = lambda path: json.loads(Path(path).read_text(encoding='utf-8'))
    summarize(load(args.plan), [load(p) for p in Path(args.predictions).glob('*.json')],
              [load(p) for p in Path(args.runtimes).glob('*/inputs.json')],
              args.check_only.lower() == 'true', args.outdir, args.rna_qc_link)


if __name__ == '__main__':
    main()
