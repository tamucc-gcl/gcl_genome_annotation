#!/usr/bin/env python3
"""Validate final identified inputs and prepare search-only protein sequences."""
import argparse
import csv
import json
from pathlib import Path

from model_metrics import sequences
from prepare_protein_evidence import checksum

TOOLS = ('interproscan', 'eggnog', 'diamond', 'kofamscan', 'deepsig')


def read_table(path, key):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    if not fields or key not in fields or any(None in r or not r.get(key) for r in rows):
        raise ValueError(f'Malformed table: {path}')
    result = {r[key]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f'Duplicate {key}: {path}')
    return fields, result


def prepare(aid, proteins, gff, assessment, validation, outdir):
    state = json.loads(Path(validation).read_text())
    if state['status'] != 'ready' or state['stage'] != 'identified_annotation' or state['assembly_id'] != aid:
        raise ValueError('Expected accepted identified annotation')
    for path, key in [(proteins, 'proteins_sha256'), (gff, 'gff3_sha256')]:
        if checksum(path) != state[key]:
            raise ValueError(f'Identified annotation changed: {path}')
    if checksum(assessment) != state['outputs']['transcript_assessment.tsv']:
        raise ValueError('Identified assessment changed')
    seqs = sequences(proteins)
    _, metrics = read_table(assessment, 'transcript_id')
    if set(seqs) != set(metrics) or len(seqs) != state['coding_transcripts']:
        raise ValueError('Protein/assessment IDs differ')
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    rows = []
    with (out/'eligible.fa').open('w') as handle:
        for tid, seq in sorted(seqs.items()):
            body = seq[:-1] if seq.endswith('*') else seq
            reason = ('empty_after_terminal_stop' if not body else
                      'internal_stop' if '*' in body else
                      'unsupported_residue' if set(body)-set('ACDEFGHIKLMNPQRSTVWYBXZJUO') else '.')
            eligible = reason == '.'
            rows.append(dict(protein_id=tid, eligible=str(eligible).lower(), reason=reason,
                             terminal_stop_removed=str(seq.endswith('*')).lower()))
            if eligible:
                handle.write(f'>{tid}\n{body}\n')
    with (out/'eligibility.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    result = dict(status='ready', assembly_id=aid, proteins=len(seqs),
                  eligible=sum(r['eligible']=='true' for r in rows),
                  source_proteins_sha256=checksum(proteins), source_gff3_sha256=checksum(gff))
    (out/'inputs.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','proteins','gff','assessment','validation','outdir'):
        parser.add_argument('--'+key, required=True)
    a = parser.parse_args()
    prepare(a.assembly_id, a.proteins, a.gff, a.assessment, a.validation, a.outdir)
