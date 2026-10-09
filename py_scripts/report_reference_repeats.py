#!/usr/bin/env python3
"""Reference/repeat status report independent of expensive scientific processing."""
import argparse
import json
from pathlib import Path


def summarize(preflight, references, libraries, repeats, outdir, benchmarks=()):
    refs, libs, masks = ({row['assembly_id']: row for row in rows} for rows in (references, libraries, repeats))
    comparisons = {row['assembly_id']: row for row in benchmarks}
    rows = []
    for request in preflight['assemblies']:
        aid = request['assembly_id']
        ref, lib, mask = refs.get(aid, {}), libs.get(aid, {}), masks.get(aid, {})
        benchmark = comparisons.get(aid)
        failed = next((state for state in (ref, lib, mask, benchmark or {}) if state.get('status') == 'failed'), None)
        ready = ref.get('status') == 'ready' and lib.get('status') == 'ready' and mask.get('status') == 'ready'
        ready = ready and (benchmark is None or benchmark.get('status') == 'ready')
        rows.append({'assembly_id': aid, 'status': 'ready' if ready else 'blocked',
                     'reason': failed.get('reason') if failed else (None if ready else 'Missing required reference/repeat result'),
                     'reference': ref, 'library': lib, 'repeats': mask, 'benchmark': benchmark})
    result = {'schema_version': 1, 'stage': 'reference_repeats',
              'status': 'ready' if all(row['status'] == 'ready' for row in rows) else 'incomplete',
              'validation_scope': 'Normalized reference, whole-contig masking, case-only identity and repeat coordinates',
              'assemblies': rows}
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    (outdir / 'reference_repeat_status.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# Reference and repeat preparation', '', f"Status: **{result['status']}**", '',
             '| Assembly | Status | Retained contigs | Retained bases | Repeat union bases | Detail |',
             '|---|---|---:|---:|---:|---|']
    for row in rows:
        aid, stats, mask = row['assembly_id'], row['reference'].get('stats', {}), row['repeats']
        lines.append(f"| {aid} | {row['status']} | {stats.get('retained_contigs', 'NA')} | "
                     f"{stats.get('retained_bases', 'NA')} | {mask.get('repeat_interval_union_bases', 'NA')} | "
                     f"[reference](../filtered_assemblies/{aid}/{aid}.reference.json) |")
        if row['reason']:
            lines += ['', f"{aid}: {row['reason'].replace('|', '/').replace(chr(10), ' ')}", '']
    lines += ['', 'Masking batches contain whole contigs, with no artificial tile boundaries.',
              'Accepted masking preserves contig order, IDs, lengths and uppercase sequence.',
              'Repeat coverage uses interval unions; class coverage can overlap and must not be summed.',
              'The legacy .tbl category lengths are alignment-fragment sums, not unique coverage.', '',
              'Normalized FASTA, its index/dictionary and excluded-contig list are in `../filtered_assemblies/`.',
              'Repeat tracks, masked FASTA and native tool summaries are in `../repeats/`.',
              '[Resolved statuses, identities and coverage](reference_repeat_status.json)', '',
              'This entry stops before RNA alignment and gene prediction.', '']
    (outdir / 'README.md').write_text('\n'.join(lines))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('preflight', 'references', 'libraries', 'repeats', 'outdir'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--benchmarks')
    args = parser.parse_args()
    rows = []
    for folder in (args.references, args.libraries, args.repeats):
        Path(folder).mkdir(exist_ok=True)
        rows.append([json.loads(path.read_text()) for path in Path(folder).glob('*.json')])
    benchmarks = [json.loads(path.read_text()) for path in Path(args.benchmarks).glob('*.json')] if args.benchmarks else []
    summarize(json.loads(Path(args.preflight).read_text()), *rows, args.outdir, benchmarks)


if __name__ == '__main__':
    main()
