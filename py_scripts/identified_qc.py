#!/usr/bin/env python3
"""Reuse baseline representative selection/reporting for the identified selected set."""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote
from annotation_baseline import prepare, report, decode
from model_metrics import sequences
from prepare_protein_evidence import checksum


def prepare_selected(aid,settings,state,gff,proteins,outdir,genome=None):
    if state.get('stage')!='identified_annotation':raise ValueError('Selected QC requires validated final identifiers')
    data=prepare(aid,settings,state,gff,proteins,outdir,genome=genome)
    # Preserve the accepted before-QC tie rule through ID changes: equal-length
    # isoforms are ranked by their original ID, not their newly numbered ID.
    seqs=sequences(proteins);tie=state['representative_tie_ids'];groups=defaultdict(list)
    if set(tie)!=set(seqs):raise ValueError('Representative tie provenance differs from final proteins')
    for line in Path(gff).read_text().splitlines():
        if not line or line.startswith('#'):continue
        c=line.split('\t')
        if c[2] in ('mRNA','transcript'):
            attrs=dict(v.split('=',1) for v in c[8].split(';') if '=' in v)
            groups[unquote(attrs['Parent'])].append(unquote(attrs['ID']))
    with Path(outdir,f'{aid}.fa').open('w') as fasta,Path(outdir,'representatives.tsv').open('w') as table:
        table.write('gene_id\ttranscript_id\tprotein_length\tisoforms\n')
        for gid,tids in sorted(groups.items()):
            tid=min(tids,key=lambda t:(-len(seqs[t].rstrip('*')),tie[t]))
            fasta.write(f'>{tid}\n{seqs[tid]}\n');table.write(f'{gid}\t{tid}\t{len(seqs[tid].rstrip("*"))}\t{len(tids)}\n')
    data.update(model_set='selected',selection='Longest protein per coding gene; ties by original transcript ID through final ID map',representative_sha256=checksum(Path(outdir,f'{aid}.fa')))
    Path(outdir,'selection.json').write_text(json.dumps(data,indent=2)+'\n')
    return data

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','report'))
    p.add_argument('--request',required=True);p.add_argument('--outdir',required=True);a=p.parse_args();r=decode(a.request)
    if a.action=='prepare':
        state=json.loads(Path(r['validation']).read_text())
        prepare_selected(r['assembly_id'],r['settings'],state,r['gff'],r['proteins'],a.outdir,genome=r['genome'])
    else:
        data=report(r['selections'],r['agat'],r['busco'],r['omark'],a.outdir)
        data.update(stage='selected_annotation_qc',scope='Selected identified models: AGAT all isoforms; BUSCO/OMArk one longest protein per gene')
        Path(a.outdir,'baseline_status.json').write_text(json.dumps(data,indent=2)+'\n')
        path=Path(a.outdir,'README.md');path.write_text(path.read_text().replace('# Unfiltered annotation QC','# Selected annotation QC').replace('Unfiltered coding models:','Selected identified coding models:').replace('original BRAKER models','selected identified models'))
