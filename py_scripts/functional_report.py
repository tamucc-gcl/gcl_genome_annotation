#!/usr/bin/env python3
"""Validate functional decoration without changing model selection or confidence."""
import argparse
import json
from pathlib import Path
from urllib.parse import unquote

from functional_inputs import read_table
from prepare_protein_evidence import checksum


def geometry(path):
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'): continue
            c = line.rstrip('\n').split('\t')
            if len(c) != 9: raise ValueError('Malformed annotated GFF')
            a = dict(x.split('=',1) for x in c[8].split(';') if '=' in x)
            rows.append(tuple(c[:8])+tuple(a.get(k,'.') for k in ('ID','Parent')))
    return sorted(rows)


def report(aid, original, annotated, table, merge_status, eligibility, outdir):
    if geometry(original) != geometry(annotated):
        raise ValueError('Functional decoration changed model geometry, phase, source or IDs')
    state = json.loads(Path(merge_status).read_text())
    if state['status'] != 'ready' or state['assembly_id'] != aid or checksum(table) != state['table_sha256']:
        raise ValueError('Functional merge provenance mismatch')
    _, rows = read_table(table, 'protein_id')
    products = {}
    for line in Path(annotated).read_text().splitlines():
        if not line or line.startswith('#'): continue
        c = line.split('\t')
        if c[2] in ('mRNA','transcript'):
            a = dict(x.split('=',1) for x in c[8].split(';') if '=' in x)
            if a['ID'] in products: raise ValueError('Duplicate decorated transcript')
            products[a['ID']] = unquote(a['product'])
    if set(products) != set(rows): raise ValueError('Decorated model/table IDs differ')
    named = sum(value != 'hypothetical protein' for value in products.values())
    state.update(stage='functional_annotation', coding_structure_unchanged=True, models_removed=0,
                 named_transcripts=named, hypothetical_transcripts=len(rows)-named,
                 original_gff3_sha256=checksum(original), annotated_gff3_sha256=checksum(annotated),
                 eligibility_sha256=checksum(eligibility),
                 scope='Functional attributes on every selected identified transcript; accepted assessment joined unchanged; no new filtering or combined confidence score')
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    (out/'functional_status.json').write_text(json.dumps(state, indent=2)+'\n')
    lines = ['# Functional annotation', '', 'Status: **ready**', '',
             f'Genes: {state["genes"]:,}. Coding transcripts: {len(rows):,}.', '',
             f'Named transcripts: {named:,}. Hypothetical proteins: {len(rows)-named:,}.', '',
             '| Tool | Hits | No hit | Ineligible | Disabled |', '|---|---:|---:|---:|---:|']
    for tool, counts in state['tools'].items():
        lines.append('| '+tool+' | '+' | '.join(str(counts.get(k,0)) for k in ('hit','no_hit','ineligible','disabled'))+' |')
    lines += ['', 'Native results:']
    suffixes = dict(interproscan='interproscan.tsv',eggnog='eggnog.annotations',diamond='diamond.tsv',
                    kofamscan='kofamscan.tsv',deepsig='deepsig.gff3')
    for tool, counts in state['tools'].items():
        if counts.get('hit',0)+counts.get('no_hit',0):
            lines.append(f'- [{tool}]({aid}.{suffixes[tool]})')
    lines += ['', f'[Annotated coding GFF3]({aid}.annotated.gff3) | [Functional and assessment table]({aid}.functional_annotations.tsv) | '
              f'[Gene products]({aid}.gene_product_summary.tsv)', '',
              '[Search eligibility](eligibility.tsv) | [Validation and provenance](functional_status.json)', '',
              f'[Identified proteins, CDS and complete ID maps](../../structural_annotation/identified/{aid}/README.md) | '
              f'[Structural QC](../../qc/structural/{aid}/README.md)', '',
              'Only a terminal stop is removed in search inputs. Internal-stop or invalid-residue sequences remain in the annotation and table, but are ineligible for searches. '
              'Partial proteins are not excluded merely for lacking a start methionine. Original FASTAs are unchanged.', '',
              'Products follow the existing priority: DIAMOND description, eggNOG name/description, significant KofamScan KO, then hypothetical protein. '
              'DeepSig Chain records alone do not indicate signal peptides. Raw tool outputs retain all reported matches.', '',
              'RNA/protein AED, junction and TE metrics are inherited from the accepted assessment. OMArk classes apply only to assessed representative proteins; '
              'other isoforms are not classified. No further models are removed. Other genome feature tracks are not included in this coding GFF3 yet.', '']
    (out/'README.md').write_text('\n'.join(lines))
    return state


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','original','annotated','table','merge-status','eligibility','outdir'):
        p.add_argument('--'+key, required=True)
    a = p.parse_args()
    report(a.assembly_id,a.original,a.annotated,a.table,a.merge_status,a.eligibility,a.outdir)
