#!/usr/bin/env python3
"""Validate native AED model/evidence inputs independently of reporting."""
import argparse
import json
from pathlib import Path

from prepare_protein_evidence import checksum


def check(gff, validation, current_status, rna, protein, output):
    validation = json.loads(Path(validation).read_text())
    current = json.loads(Path(current_status).read_text())
    if validation.get('status') != 'ready' or current.get('status') != 'ready':
        raise ValueError('Benchmark requires accepted predictions and completed current AED')
    if checksum(gff) != validation['gff3_sha256'] or checksum(gff) != current['model_sha256'] or validation['reference_identity'] != current['reference_identity']:
        raise ValueError('Benchmark models differ from accepted prediction/current AED')
    for kind, path in [('rna',rna),('protein',protein)]:
        source = current['evidence'][kind]
        if source['available'] != bool(path):
            raise ValueError(f'{kind}: supplied evidence differs from current AED availability')
        if path and checksum(path) != source['source_sha256']:
            raise ValueError(f'{kind}: evidence changed since current AED')
    data = dict(status='ready', model_sha256=checksum(gff), reference_identity=current['reference_identity'],
                evidence=current['evidence'])
    Path(output).write_text(json.dumps(data,indent=2)+'\n')
    return data



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('gff','validation','current-status','out'):
        parser.add_argument('--'+arg,required=True)
    parser.add_argument('--rna')
    parser.add_argument('--protein')
    args = parser.parse_args()
    check(args.gff,args.validation,args.current_status,args.rna,args.protein,args.out)


if __name__=='__main__': main()
