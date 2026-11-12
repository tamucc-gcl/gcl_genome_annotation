#!/usr/bin/env python3
"""Validate native AED model/evidence inputs independently of reporting."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from urllib.parse import unquote
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


def split_proteins(source, outputs, source_roles, expected_sha, mapping):
    """Miniprot parent + child projection, retaining complete alignments by source."""
    parents = {}
    count = Counter()
    sha = hashlib.sha256()
    with Path(source).open('rb') as src:
        for raw in src:
            sha.update(raw)
            if raw.startswith(b'#') or not raw.strip():
                continue
            c = raw.decode().rstrip('\r\n').split('\t')
            if len(c) != 9:
                raise ValueError('Unexpected miniprot record')
            a = dict(x.split('=',1) for x in c[8].split(';') if '=' in x)
            if 'Target' in a:
                pid = unquote(a['Target']).split()[0]
                source_id = pid.split('__')[0]
                role = source_roles.get(source_id)
                if role is None:
                    raise ValueError('Unclassified protein Target: '+pid)
                if 'ID' not in a and 'Parent' not in a:
                    raise ValueError('Target alignment lacks ID')
                if 'ID' in a and a['ID'] in parents and parents[a['ID']] != role:
                    raise ValueError('Conflicting miniprot parent')
                if 'ID' in a:
                    parents[a['ID']] = role
                if 'Parent' in a and any(parents.get(p) != role for p in a['Parent'].split(',')):
                    raise ValueError('Child Target differs from parent source')
            elif 'Parent' in a:
                roles = {parents.get(p) for p in a['Parent'].split(',')}
                if len(roles) != 1 or None in roles:
                    raise ValueError('Unresolved/cross-source miniprot parent')
                role = roles.pop()
            else:
                raise ValueError('Alignment cannot be assigned to evidence source')
            outputs[role].write(raw)
            count[role] += 1
            if c[2] == 'mRNA' and 'Target' in a:
                mapping.writerow([a['ID'], pid, source_id, role])
    if sha.hexdigest() != expected_sha:
        raise ValueError('Existing protein alignment differs from accepted AED source')
    return dict(count)


def partition(data, protein, provenance, outdir):
    roles = {}
    if protein:
        manifest = json.loads(Path(provenance).read_text())
        for entry in manifest['sources']:
            kind = entry['request']['kind']
            if kind not in ('orthodb','ncbi','local'):
                raise ValueError('Unsupported protein evidence kind: '+kind)
            roles[entry['source_id']] = 'orthodb' if kind == 'orthodb' else 'reference'
        if not roles:
            raise ValueError('Protein alignment lacks source provenance')
    outdir = Path(outdir)
    with (outdir/'reference.gff').open('wb') as ref, (outdir/'orthodb.gff').open('wb') as db, (outdir/'protein_ids.tsv').open('w',newline='') as ids:
        writer = csv.writer(ids, delimiter='\t')
        writer.writerow(['alignment_id','protein_id','source_id','role'])
        if protein:
            split_proteins(protein, {'reference':ref,'orthodb':db}, roles, data['evidence']['protein']['source_sha256'], writer)
    for role, filename in [('reference','inputs.json'),('orthodb','orthodb.json')]:
        result = dict(data, source_roles=roles)
        result['evidence'] = dict(data['evidence'], protein=dict(data['evidence']['protein'], available=role in roles.values(), pool=role))
        (outdir/filename).write_text(json.dumps(result, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('gff','validation','out'):
        parser.add_argument('--'+arg,required=True)
    parser.add_argument('--rna')
    parser.add_argument('--protein')
    parser.add_argument('--protein-provenance')
    parser.add_argument('--split-protein-pools', action='store_true')
    args = parser.parse_args()
    if args.split_protein_pools and args.protein and not args.protein_provenance:
        parser.error('--protein-provenance is required with --protein')
    data = check(args.gff,args.validation,args.rna,args.protein,args.out)
    if args.split_protein_pools:
        partition(data,args.protein,args.protein_provenance,Path(args.out).parent)


if __name__=='__main__': main()
