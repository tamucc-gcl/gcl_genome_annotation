#!/usr/bin/env python3
"""Validate native AED model/evidence inputs independently of reporting."""
import argparse
import json
from pathlib import Path

from prepare_protein_evidence import checksum


def check(gff, validation, rna, protein, output):
    validation = json.loads(Path(validation).read_text())
    if validation.get('status') != 'ready' or validation.get('stage') != 'final_annotation':
        raise ValueError('Native AED requires a validated final annotation')
    mode = validation['mode']
    if mode not in ('rna_protein', 'rna', 'protein', 'ab_initio'):
        raise ValueError('Unknown final annotation evidence mode')
    if checksum(gff) != validation['gff3_sha256']:
        raise ValueError('Native AED models differ from accepted annotation')
    available = dict(rna=mode in ('rna', 'rna_protein'), protein=mode in ('protein', 'rna_protein'))
    evidence = {}
    for kind, path in [('rna', rna), ('protein', protein)]:
        if available[kind] != bool(path):
            raise ValueError(f'{kind}: supplied evidence differs from annotation mode availability')
        evidence[kind] = dict(available=available[kind], source_sha256=checksum(path) if path else None)
    data = dict(status='ready', model_sha256=checksum(gff),
                reference_identity=validation['reference_identity'], evidence=evidence)
    Path(output).write_text(json.dumps(data, indent=2)+'\n')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('gff','validation','out'):
        parser.add_argument('--'+arg,required=True)
    parser.add_argument('--rna')
    parser.add_argument('--protein')
    args = parser.parse_args()
    check(args.gff,args.validation,args.rna,args.protein,args.out)


if __name__=='__main__': main()
