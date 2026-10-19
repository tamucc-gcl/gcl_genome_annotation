#!/usr/bin/env python3
"""Native TE database/search execution, independent of flag/report settings."""
import argparse
import json
from pathlib import Path
import subprocess

from prepare_protein_evidence import checksum, records


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def database(library,translated,folder):
    with Path(translated).open() as handle:seqs=list(records(handle))
    ids=[h.split()[0] for h,s in seqs]
    if len(ids)!=len(set(ids)) or any(not s for h,s in seqs):raise ValueError('Duplicate/empty translated repeat ID')
    result=dict(status='ready',library_sha256=checksum(library),translated_sha256=checksum(translated),
        translated_sequences=len(ids),assessed=bool(ids),translation='SeqKit six frames; stop symbols converted to X; minimum 30 aa; standard code',
        reason=None if ids else 'Empty repeat library or no translated sequence >=30 aa')
    if ids and not (Path(folder)/'repeats.dmnd').is_file():raise ValueError('Missing translated repeat database')
    dump(Path(folder)/'database.json',result)


def search(proteins,db,threads,evalue,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True);db=Path(db)
    source=json.loads((db/'database.json').read_text())
    if source['status']!='ready':raise ValueError('TE database not ready')
    hits=out/'hits.tsv'
    if source['assessed']:
        subprocess.run(['diamond','blastp','--db',str(db/'repeats.dmnd'),'--query',str(proteins),'--out',str(hits),
            '--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','evalue','bitscore','qstart','qend',
            '--evalue',str(evalue),'--max-target-seqs','1','--threads',str(threads)],check=True)
        if not hits.is_file():raise ValueError('DIAMOND succeeded without an output file')
    else:hits.write_text('')
    dump(out/'search.json',dict(status='ready',assessed=source['assessed'],reason=source['reason'],
        proteins_sha256=checksum(proteins),library_sha256=source['library_sha256'],evalue=evalue,
        hits_sha256=checksum(hits),hits=sum(bool(line.strip()) for line in hits.open()),database=source))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('database')
    for key in ('library','translated','folder'):a.add_argument('--'+key,required=True)
    a=sub.add_parser('search')
    for key in ('proteins','db','outdir'):a.add_argument('--'+key,required=True)
    a.add_argument('--threads',type=int,required=True);a.add_argument('--evalue',type=float,required=True)
    a=p.parse_args()
    if a.action=='database':database(a.library,a.translated,a.folder)
    else:search(a.proteins,a.db,a.threads,a.evalue,a.outdir)
