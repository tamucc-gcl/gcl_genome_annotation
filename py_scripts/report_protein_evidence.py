#!/usr/bin/env python3
"""Report assembled prediction protein evidence without changing source tasks."""
import argparse
import base64
import json
import re
from pathlib import Path
from urllib.parse import urlsplit

SCHEMA = 1


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def protein_source_link(source):
    request = source.get('request', {})
    kind = request.get('kind')
    if kind == 'orthodb':
        clade = request.get('clade')
        label = ' '.join(part for part in ('OrthoDB', clade) if part)
        if request.get('release'):
            label += f" ({request['release']})"
        archive = urlsplit(request.get('url') or '')
        # The release catalogue lists the partition alongside its provider checksum.
        directory, _, filename = archive.path.rpartition('/')
        if (archive.scheme in ('http', 'https') and archive.netloc
                and clade and filename == f'{clade}.fa.gz'):
            page = archive._replace(path=directory + '/', query='', fragment='').geturl()
            return f'[{label}]({page})'
        return label
    if kind == 'ncbi':
        accession = request.get('accession')
        if accession and re.fullmatch(r'GC[AF]_[0-9]+\.[0-9]+', accession):
            return f'[NCBI {accession}](https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/)'
        return f'NCBI {accession or "proteome"}'
    if kind == 'local':
        path = request.get('local_path')
        return f"Local FASTA {Path(path).name}" if path else 'Local FASTA'
    return kind or source['source_id']


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
        for source in sorted(sources.values(), key=lambda source: source['source_id']):
            lines.append(f"| {escape(protein_source_link(source))} | "
                         f"{escape(source['request']['kind'])} | {escape(source['status'])} |")
    lines += ['', '[Resolved protein evidence](protein_status.json)', '',
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
