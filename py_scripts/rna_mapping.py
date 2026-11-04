#!/usr/bin/env python3
"""Stable mapping-only contracts; strand policy changes must not invalidate mapping."""
import argparse
import base64
import hashlib
import json
from pathlib import Path


def selected_reads(spec, paths):
    from prepare_rna import digest
    expected = spec.get('sha256') or []
    if len(paths) != (2 if spec['layout'] == 'paired' else 1) or len(expected) != len(paths):
        raise ValueError('Selected-read checksum/file count mismatch')
    if [digest(path) for path in paths] != expected:
        raise ValueError('FASTQ bytes differ from the accepted selected-read manifest; rerun RNA preparation')


def reference(source, identity):
    # Reuse the accepted reference parser and content identity contract.
    from reference_repeats import fasta, sha
    rows = [{'id': header.split()[0], 'length': len(seq),
             'uppercase_sha256': hashlib.sha256(seq.upper().encode('ascii')).hexdigest()}
            for header, seq in fasta(source)]
    if rows != identity['contigs']:
        raise ValueError('Indexed reference differs from the accepted contig identities')
    value = hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if value != identity['reference_identity']:
        raise ValueError('Invalid reference identity checksum')
    return {'contigs': rows, 'reference_identity': value, 'indexed_fasta_sha256': sha(source)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['reads', 'reference'])
    for name in ('request', 'reference', 'identity', 'output'):
        parser.add_argument('--' + name)
    args = parser.parse_args()
    if args.action == 'reads':
        request = json.loads(base64.b64decode(args.request))
        selected_reads(request['spec'], request['paths'])
    else:
        result = reference(args.reference, json.loads(Path(args.identity).read_text(encoding='utf-8')))
        Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
