#!/usr/bin/env python3
"""Adapt FLAG_TE_GENES overlap/homology evidence; retain every candidate model."""
import argparse
import csv
import json
import math
from pathlib import Path
from urllib.parse import unquote, quote
from collections import defaultdict

from prepare_protein_evidence import checksum


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def repeat_bed(source,identity,output):
    ref=json.loads(Path(identity).read_text());lengths={c['id']:c['length'] for c in ref['contigs']}
    with Path(source).open() as handle,Path(output).open('w') as out:
        for line in handle:
            if not line.strip() or line.startswith('#'):continue
            r=line.rstrip('\n').split('\t')
            if len(r)!=9:raise ValueError('Repeat annotation must have nine GFF columns')
            s,e=int(r[3])-1,int(r[4])
            if r[0] not in lengths or not 0<=s<e<=lengths[r[0]]:raise ValueError('Repeat outside accepted reference')
            out.write(f'{r[0]}\t{s}\t{e}\n')


def flag(aid,gff,metrics,coverage,search_dir,overlap,outdir):
    if not 0<overlap<=1:raise ValueError('TE overlap threshold must be >0 and <=1')
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    data=json.loads(Path(metrics).read_text());search_dir=Path(search_dir);native=json.loads((search_dir/'search.json').read_text())
    if data['status']!='ready' or checksum(gff)!=data['model_sha256']:raise ValueError('TE models differ from metrics')
    if native['status']!='ready' or native['proteins_sha256']!=data['proteins_sha256'] or checksum(search_dir/'hits.tsv')!=native['hits_sha256']:
        raise ValueError('TE protein search provenance mismatch')
    tx={r['transcript_id']:r for r in data['metrics']};genes={r['gene_id'] for r in tx.values()};totals=defaultdict(lambda:[0,0]);seen=set()
    with Path(coverage).open() as handle:
        for line in handle:
            r=line.rstrip('\n').split('\t')
            if len(r)!=10:raise ValueError('Unexpected BEDTools coverage columns')
            kind,gid=r[3].split('|',1);gid=unquote(gid);start,end=int(r[1]),int(r[2]);bp,total=int(r[7]),int(r[8])
            if gid not in genes or kind not in ('gene','cds','exon') or not 0<=bp<=total or total!=end-start:
                raise ValueError('Invalid gene repeat coverage')
            key=(gid,kind,r[0],start,end)
            if key in seen:raise ValueError('Duplicate gene coverage block')
            seen.add(key);totals[gid,kind][0]+=bp;totals[gid,kind][1]+=total
    if set(totals)!={(gid,kind) for gid in genes for kind in ('gene','cds','exon')}:
        raise ValueError('Incomplete gene repeat coverage')
    hits={}
    with (search_dir/'hits.tsv').open() as handle:
        for line in handle:
            if not line.strip():continue
            r=line.rstrip('\n').split('\t')
            if len(r)!=10 or r[0] not in tx:raise ValueError('Unknown TE query or malformed DIAMOND hit')
            pid,target=r[:2];identity,length,qlen,slen,evalue,bits,qstart,qend=map(float,r[2:])
            if (not all(math.isfinite(v) for v in (identity,length,qlen,slen,evalue,bits,qstart,qend))
                or not 0<=identity<=100 or evalue<0 or evalue>native['evalue'] or not 1<=qstart<=qend<=qlen):
                raise ValueError('Invalid TE alignment metrics')
            hit=dict(transcript_id=pid,hit=target,pident=identity,qcov=100*(qend-qstart+1)/qlen,evalue=evalue,bitscore=bits)
            if pid not in hits or bits>hits[pid]['bitscore']:hits[pid]=hit
    by_gene=defaultdict(list)
    for pid,hit in hits.items():by_gene[tx[pid]['gene_id']].append(hit)
    rows=[];flagged=set()
    for gid in sorted(genes):
        candidates=by_gene[gid]
        best=max(candidates,key=lambda r:(r['bitscore'],r['transcript_id'])) if candidates else None
        row=dict(gene_id=gid)
        for kind in ('gene','exon','cds'):
            bp,total=totals[gid,kind];row[kind+'_repeat_bp']=bp;row[kind+'_bp']=total;row[kind+'_repeat_fraction']=bp/total
        evidence=[]
        if row['gene_repeat_fraction']>=overlap:evidence.append('repeat_overlap')
        if best:evidence.append('repeat_library_homology')
        row.update(te_flagged=bool(evidence),te_evidence=';'.join(evidence) or '.',homology_assessed=native['assessed'],
            homology_transcripts=','.join(sorted(h['transcript_id'] for h in candidates)) or '.',
            best_hit=best['hit'] if best else '.',best_pident=best['pident'] if best else None,
            best_qcov=best['qcov'] if best else None,best_evalue=best['evalue'] if best else None)
        if evidence:flagged.add(gid)
        rows.append(row)
    with (out/'te_gene_report.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({k:'.' if v is None else v for k,v in row.items()})
    with Path(gff).open() as source,(out/f'{aid}.te_flagged.gff3').open('w') as dest:
        for line in source:
            if line.startswith('#') or not line.strip():dest.write(line);continue
            r=line.rstrip('\n').split('\t');attrs=dict(x.split('=',1) for x in r[8].split(';') if '=' in x)
            gid=unquote(attrs.get('ID','')) if r[2]=='gene' else None
            owners={tx[t]['gene_id'] for t in unquote(attrs.get('Parent','')).split(',') if t in tx}
            if r[2] in ('mRNA','transcript'):owners={unquote(attrs.get('Parent',''))}
            text='putative_TE_derived' if gid in flagged else 'parent_flagged_TE_derived' if owners&flagged else None
            if text:
                # Preserve IDs/Parents/coordinates and escape new attribute values.
                note=attrs.get('Note');items=[x for x in r[8].split(';') if x and not x.startswith('Note=')]
                items.append('Note='+(note+',' if note else '')+quote(text,safe=''))
                r[8]=';'.join(items);dest.write('\t'.join(r)+'\n')
            else:dest.write(line)
    status=dict(status='ready',stage='te_flags',assembly_id=aid,model_sha256=data['model_sha256'],
        reference_identity=data['reference_identity'],genes=len(genes),flagged_genes=len(flagged),overlap_threshold=overlap,
        overlap_scope='Union of repeat intervals over gene span; exon/CDS gene-union fractions also reported',
        search=native,transcript_hits=hits,report_sha256=checksum(out/'te_gene_report.tsv'),
        flagged_gff3_sha256=checksum(out/f'{aid}.te_flagged.gff3'),scope='Diagnostic repeat-derived flags; no model removal; homology to pre-filter repeat library can include host-gene sequence')
    dump(out/'te_status.json',status);return status


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('repeat-bed')
    for key in ('source','identity','out'):a.add_argument('--'+key,required=True)
    a=sub.add_parser('flag')
    for key in ('assembly-id','gff','metrics','coverage','search-dir','outdir'):a.add_argument('--'+key,required=True)
    a.add_argument('--overlap',type=float,required=True)
    a=p.parse_args()
    if a.action=='repeat-bed':repeat_bed(a.source,a.identity,a.out)
    else:flag(a.assembly_id,a.gff,a.metrics,a.coverage,a.search_dir,a.overlap,a.outdir)
