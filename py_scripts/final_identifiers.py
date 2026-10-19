#!/usr/bin/env python3
"""Adapt RENAME_GENE_MODELS' assembly_gN/.tN scheme with complete, validated maps."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
from urllib.parse import quote, unquote

from evidence_aed import models
from model_metrics import sequences
from model_selection import table, write_table
from prepare_protein_evidence import checksum


def rename(aid,selection,identity,outdir,prefix=None):
    source=Path(selection);out=Path(outdir)
    status=json.loads((source/'selection_status.json').read_text())
    ref=json.loads(Path(identity).read_text());lengths={r['id']:r['length'] for r in ref['contigs']}
    if (status.get('stage')!='model_selection' or status.get('status')!='ready' or status.get('assembly_id')!=aid or
        status['reference_identity']!=ref['reference_identity']):raise ValueError('Identifiers require accepted selection/reference')
    for name,sha in status['outputs'].items():
        if checksum(source/name)!=sha:raise ValueError('Selection output changed: '+name)
    if not status['retained_transcripts']:raise ValueError('No retained coding models; review selection before final identifiers/QC')
    prefix=prefix or aid
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*',prefix):raise ValueError('Final ID prefix must use letters, digits, underscore, dot or hyphen')
    gff=source/f'retained/{aid}.gff3';tx=models(gff,lengths)
    def parse(path):
        result=[]
        for ordinal,line in enumerate(Path(path).read_text().splitlines(),1):
            if not line or line.startswith('#'):continue
            c=line.split('\t')
            if len(c)!=9:raise ValueError('Malformed selected GFF')
            a=dict(v.split('=',1) for v in c[8].split(';') if '=' in v)
            result.append((ordinal,c,a,unquote(a.get('ID','')),tuple(unquote(p) for p in a.get('Parent','').split(',') if p)))
        return result
    rows=parse(gff);nodes={};order={r['id']:i for i,r in enumerate(ref['contigs'])}
    for ordinal,c,a,fid,parents in rows:
        if fid:
            if fid in nodes and (nodes[fid][1][2],nodes[fid][4])!=(c[2],parents):raise ValueError('Inconsistent multipart ID')
            nodes.setdefault(fid,(ordinal,c,a,fid,parents))
    key=lambda r:(order[r[1][0]],int(r[1][3]),int(r[1][4]),r[1][6],r[3])
    genes=sorted((r for r in nodes.values() if r[1][2]=='gene'),key=key)
    mapping={r[3]:f'{prefix}_g{i}' for i,r in enumerate(genes,1)};gmap=dict(mapping)
    tmap={};counts=defaultdict(int)
    for row in sorted((r for r in nodes.values() if r[1][2] in ('mRNA','transcript')),key=lambda r:(gmap[r[4][0]],key(r))):
        old=row[3];parents=row[4]
        if len(parents)!=1 or parents[0] not in gmap:raise ValueError('Unresolved transcript parent')
        parent=gmap[parents[0]];counts[parent]+=1;tmap[old]=f'{parent}.t{counts[parent]}'
    mapping.update(tmap);child_counts=defaultdict(int)
    pending=[r for r in rows if r[1][2] not in ('gene','mRNA','transcript')];anonymous={}
    while pending:
        remaining=[];progress=False
        for ordinal,c,a,fid,parents in sorted(pending,key=key):
            if not parents or any(p not in mapping for p in parents):remaining.append((ordinal,c,a,fid,parents));continue
            if fid and fid in mapping:progress=True;continue # multipart feature retains one ID
            parent=sorted(mapping[p] for p in parents)[0];kind=re.sub(r'[^A-Za-z0-9_]','_',c[2])
            child_counts[(parent,kind)]+=1;new=f'{parent}.{kind}{child_counts[(parent,kind)]}'
            if fid:mapping[fid]=new
            else:anonymous[ordinal]=new
            progress=True
        if not progress:raise ValueError('Unresolved/cyclic child parents during renaming')
        pending=remaining
    if len(set(mapping.values())|set(anonymous.values()))!=len(mapping)+len(anonymous):raise ValueError('Final IDs collide')
    out.mkdir(parents=True,exist_ok=True);map_rows=[]
    with (out/f'{aid}.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        # Parent-before-child output, preserving original feature geometry/phases.
        emitted=set();pending=list(rows)
        while pending:
            remaining=[];progress=False
            for ordinal,c,a,fid,parents in pending:
                if any(p not in emitted for p in parents):remaining.append((ordinal,c,a,fid,parents));continue
                new=mapping[fid] if fid else anonymous[ordinal]
                attrs=dict(a);attrs['ID']=quote(new,safe='._:-')
                if parents:attrs['Parent']=','.join(quote(mapping[p],safe='._:-') for p in parents)
                if 'Derives_from' in attrs:
                    links=[unquote(p) for p in attrs['Derives_from'].split(',')]
                    if any(p not in mapping for p in links):raise ValueError('Unresolved Derives_from')
                    attrs['Derives_from']=','.join(quote(mapping[p],safe='._:-') for p in links)
                for attr in ('Name','gene_id','transcript_id'):
                    if attr in attrs and unquote(attrs[attr]) in mapping:attrs[attr]=quote(mapping[unquote(attrs[attr])],safe='._:-')
                handle.write('\t'.join(c[:8]+[';'.join(k+'='+v for k,v in attrs.items())])+'\n')
                map_rows.append(dict(partition='retained',source_line=ordinal,feature_type=c[2],old_id=fid or '.',new_id=new,
                    old_parent=','.join(parents) or '.',new_parent=','.join(mapping[p] for p in parents) or '.',
                    seqid=c[0],start=c[3],end=c[4],strand=c[6],phase=c[7]))
                if fid:emitted.add(fid)
                progress=True
            if not progress:raise ValueError('Cannot serialize parent hierarchy')
            pending=remaining
    after=models(out/f'{aid}.gff3',lengths)
    if set(after)!=set(tmap.values()):raise ValueError('Renamed transcript set differs')
    for tid,m in tx.items():
        n=after[tmap[tid]]
        if any(n[k]!=m[k] for k in ('seq','start','end','strand','exons','cds')) or n['gene_id']!=gmap[m['gene_id']]:raise ValueError('Renaming changed model geometry')
    for suffix in ('proteins','cds'):
        seqs=sequences(source/f'retained/{aid}_{suffix}.fa')
        if set(seqs)!=set(tmap):raise ValueError('Selected FASTA IDs differ')
        with (out/f'{aid}_{suffix}.fa').open('w') as handle:
            for tid in sorted(seqs,key=lambda t:tmap[t]):handle.write(f'>{tmap[tid]}\n{seqs[tid]}\n')
        if sequences(out/f'{aid}_{suffix}.fa')!={tmap[t]:v for t,v in seqs.items()}:raise ValueError('Coding sequences changed')
    fields,assessment=table(source/'retained/transcript_assessment.tsv','transcript_id')
    for row in assessment:
        row.update(original_gene_id=row['gene_id'],original_transcript_id=row['transcript_id'])
        row['gene_id']=gmap[row['gene_id']];row['transcript_id']=tmap[row['transcript_id']]
    write_table(out/'transcript_assessment.tsv',fields+['original_gene_id','original_transcript_id'],assessment)
    for label,key,mp in [('transcript','transcript_id',tmap),('gene','gene_id',gmap)]:
        fields,decisions=table(source/f'{label}_decisions.tsv',key)
        for row in decisions:
            row['final_gene_id']=gmap.get(row['gene_id'],'.')
            if label=='transcript':row['final_transcript_id']=mp.get(row[key],'.')
        write_table(out/f'{label}_id_map.tsv',fields+['final_gene_id']+(['final_transcript_id'] if label=='transcript' else []),decisions)
    for ordinal,c,a,fid,parents in parse(source/f'excluded/{aid}.gff3'):
        map_rows.append(dict(partition='excluded',source_line=ordinal,feature_type=c[2],old_id=fid or '.',new_id=mapping.get(fid,'.'),
            old_parent=','.join(parents) or '.',new_parent=','.join(mapping[p] for p in parents if p in mapping) or '.',
            seqid=c[0],start=c[3],end=c[4],strand=c[6],phase=c[7]))
    write_table(out/'feature_id_map.tsv',list(map_rows[0]),map_rows)
    result=dict(status='ready',stage='identified_annotation',model_set='selected',assembly_id=aid,reference_identity=ref['reference_identity'],
        genes=len(gmap),coding_transcripts=len(tmap),prefix=prefix,selection_sha256=checksum(source/'selection_status.json'),
        original_model_sha256=checksum(gff),gff3_sha256=checksum(out/f'{aid}.gff3'),proteins_sha256=checksum(out/f'{aid}_proteins.fa'),cds_sha256=checksum(out/f'{aid}_cds.fa'),
        coordinates_and_phases_unchanged=True,coding_sequences_unchanged=True,representative_tie_ids={new:old for old,new in tmap.items()},
        scope='Deterministic reference/coordinate ordered assembly_gN/.tN IDs within this selected model set; excluded IDs map to unavailable',
        outputs={p.name:checksum(p) for p in out.iterdir() if p.is_file() and p.name not in ('annotation_status.json','README.md')})
    (out/'annotation_status.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'README.md').write_text('# Identified selected annotation\n\nStatus: **ready**\n\n'
        f"Genes: {len(gmap)}. Coding transcripts: {len(tmap)}. Prefix: `{prefix}`.\n\n"
        f'[GFF3]({aid}.gff3) | [Proteins]({aid}_proteins.fa) | [CDS]({aid}_cds.fa) | [Assessment](transcript_assessment.tsv)\n\n'
        '[All gene decisions and ID mappings](gene_id_map.tsv) | [All transcript decisions and ID mappings](transcript_id_map.tsv) | '
        '[Every retained/excluded feature record](feature_id_map.tsv) | [Validation](annotation_status.json)\n\n'
        'Removed features have no final ID unless shared with a retained feature. External evidence IDs are preserved. '
        'Coordinates, phases and coding sequences are unchanged. Numbering is deterministic for a fixed selected set; '
        'changing that set can change numbering. Keep these maps with the annotation.\n\n'
        f'[Structural QC report](../../../qc/structural/{aid}/README.md)\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','selection','identity','outdir'):p.add_argument('--'+key,required=True)
    p.add_argument('--prefix');a=p.parse_args();rename(a.assembly_id,a.selection,a.identity,a.outdir,a.prefix)
