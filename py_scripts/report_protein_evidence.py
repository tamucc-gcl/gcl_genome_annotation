#!/usr/bin/env python3
"""Report assembled prediction protein evidence without changing source tasks."""
import argparse
import base64
import json
from pathlib import Path

SCHEMA = 1


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def summarize(directories, status_path, report_path):
    assemblies = sorted((read_json(Path(folder) / 'evidence.json') for folder in directories),
                        key=lambda x: x['assembly_id'])
    result = {'schema_version': SCHEMA, 'stage': 'protein_preparation',
              'status': 'ready' if all(x['status'] in ('ready', 'not_requested') for x in assemblies) else 'incomplete',
              'assemblies': assemblies}
    write_json(status_path, result)
    def escape(value):
        return str(value).replace('|', '\\|').replace('\n', ' ').replace('\r', ' ')
    lines = ['# Prediction protein evidence', '', f"Status: **{result['status']}**", '',
             '| Assembly | Status | Proteins | Reason |', '|---|---|---|---|']
    for row in assemblies:
        label = f"[{escape(row['assembly_id'])}](assemblies/{row['assembly_id']}/evidence.json)"
        lines.append('| ' + label + ' | ' + ' | '.join(escape(row.get(key, '')) for key in
            ('status', 'protein_count', 'reason')) + ' |')
    sources = {source['source_id']: source for row in assemblies for source in row['sources']}
    if sources:
        lines += ['', '| Source | Kind | Status |', '|---|---|---|']
        for sid, source in sorted(sources.items()):
            lines.append(f"| [{sid}](../protein_sources/{sid}/source.json) | "
                         f"{escape(source['request']['kind'])} | {escape(source['status'])} |")
    lines += ['', '[Resolved protein evidence](protein_status.json)', '',
              'Source bundles: ../protein_sources/. Per-assembly outputs: assemblies/<assembly_id>/.',
              'No prediction or RNA processing was run by this entry point.', '']
    Path(report_path).write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundles-base64', required=True)
    parser.add_argument('--status', required=True)
    parser.add_argument('--report', required=True)
    args = parser.parse_args()
    directories = json.loads(base64.b64decode(args.bundles_base64).decode())
    summarize(directories, args.status, args.report)


if __name__ == '__main__':
    main()
