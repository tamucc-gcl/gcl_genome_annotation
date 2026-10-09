#!/usr/bin/env python3
"""Portable RNA preparation report; report changes do not invalidate read processing."""
import argparse
import base64
import json
import subprocess
from pathlib import Path


def summarize(plan, results, outdir, read_paths=None):
    read_paths = read_paths or {}
    states = {row['readset_id']: row for row in results}
    rows = []
    for request in plan['readsets']:
        rid = request['readset_id']
        result = states.get(rid, {'readset_id': rid, 'status': 'failed', 'reason': 'Missing preparation result'})
        rows.append({**result, 'library_id': request['library_id'], 'sample_id': request.get('sample_id'),
                     'strandedness': request['strandedness'],
                     'strand_status': 'pending_reference_alignment' if request['strandedness'] == 'auto' else 'supplied',
                     'read_paths': read_paths.get(rid, [])
                                  if result['status'] == 'ready' else []})
    status = 'ready' if all(row['status'] == 'ready' for row in rows) else 'incomplete'
    assemblies = []
    for assembly in plan['assemblies']:
        blocked = [rid for rid in assembly['rna_readsets'] if states.get(rid, {}).get('status') != 'ready']
        assemblies.append({**assembly, 'status': 'blocked' if blocked else
                           ('ready' if assembly['rna_readsets'] else 'not_requested'), 'blocked_readsets': blocked})
    result = {'schema_version': 1, 'stage': 'rna_preparation', 'status': status,
              'validation_scope': 'FASTQ integrity, preprocessing and QC; no alignment/strand inference',
              'readsets': rows, 'assemblies': assemblies}
    outdir = Path(outdir)
    (outdir / 'rna_status.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# RNA preparation', '', f'Status: **{status}**', '',
             'Read units are pairs for paired-end data and individual reads for single-end data.',
             'Automatic strandedness awaits reference-specific alignment; readiness here does not resolve it.', '',
             '| Readset | Library | Status | Selected units | Strandedness | Detail |',
             '|---|---|---|---:|---|---|']
    for row in rows:
        rid = row['readset_id']
        directory = {'raw': 'raw', 'processed': 'processed', 'selected': 'selected'}.get(row.get('stage'), 'selected')
        lines.append(f"| {rid} | {row['library_id']} | {row['status']} | {row.get('read_units', 'NA')} | "
                     f"{row['strandedness']} | [manifest]({directory}/{rid}/rna.json) |")
    for row in rows:
        if row.get('reason'):
            lines += ['', f"{row['readset_id']}: {row['reason'].replace('|', '/').replace(chr(10), ' ')}", '']
    lines += ['', '| Assembly | Status | Requested readsets |', '|---|---|---|']
    for row in assemblies:
        lines.append(f"| {row['assembly_id']} | {row['status']} | {', '.join(row['rna_readsets']) or 'none'} |")
    lines += ['', '[Raw/processed QC and trimming report](multiqc/multiqc_report.html)',
              '[Resolved settings and read identities](rna_status.json)', '',
              'Actual alignment inputs are listed in `rna_status.json`. Local and passthrough reads reuse their existing paths.',
              'Downloaded reads and newly produced cleaned/subset reads are published in the corresponding stage directory.',
              'Detailed read retention and tool commands are recorded in each processed readset manifest.', '']
    (outdir / 'README.md').write_text('\n'.join(lines))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'results', 'qc', 'fastp', 'outdir'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--read-paths', required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)
    for directory in (args.results, args.qc, args.fastp):
        Path(directory).mkdir(exist_ok=True)
    plan = json.loads(Path(args.plan).read_text())
    results = [json.loads(path.read_text()) for path in Path(args.results).glob('*.json')]
    summarize(plan, results, outdir, json.loads(base64.b64decode(args.read_paths)))
    qc = outdir / 'multiqc'
    qc.mkdir()
    if any(Path(args.qc).glob('*')) or any(Path(args.fastp).glob('*')):
        command = ['multiqc', args.qc, args.fastp, '--force', '--outdir', str(qc),
                   '--filename', 'multiqc_report.html']
        with (qc / 'multiqc.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        (qc / 'multiqc_report.html').write_text('<!doctype html><title>RNA QC</title><p>No successful RNA inputs to report.</p>')


if __name__ == '__main__':
    main()
