#!/usr/bin/env python3
"""Post-search genome-track QC; native tools retain responsibility for evidence scoring."""
import argparse
import base64
from collections import Counter
import csv
import json
import math
from pathlib import Path

from genome_tracks import dictionary
from reference_repeats import dump, sha


def read_gff(path, lengths):
    records=[]; ids=set()
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):continue
        c=line.split('\t')
        if len(c)!=9:raise ValueError('Expected nine GFF columns')
        a=dict(t.split('=',1) for t in c[8].strip(';').split(';') if '=' in t)
        if c[0] not in lengths or not 1<=int(c[3])<=int(c[4])<=lengths[c[0]]:raise ValueError('Invalid track reference coordinates')
        if not a.get('ID') or a['ID'] in ids:raise ValueError('Missing/duplicate track ID')
        ids.add(a['ID']);records.append((c,a))
    if any(p not in ids for c,a in records for p in a.get('Parent','').split(',') if p):raise ValueError('Missing track parent')
    return records


def write_gff(path, records):
    Path(path).write_text('##gff-version 3\n'+''.join('\t'.join(c[:8]+[';'.join(f'{k}={v}' for k,v in a.items())])+'\n' for c,a in records))


def trna_categories(path):
    categories={};note_index=None
    for line in Path(path).read_text().splitlines():
        cols=[s.strip() for s in line.split('\t')]
        if 'Note' in cols:note_index=cols.index('Note')
        if len(cols)<4 or not cols[1].isdigit() or not cols[2].isdigit() or not cols[3].isdigit():continue
        if note_index is None or len(cols)<=note_index:raise ValueError('Native high-confidence output has no Note column')
        start,end=int(cols[2]),int(cols[3]);key=(cols[0],min(start,end),max(start,end),'+' if start<end else '-')
        if key in categories:raise ValueError('Duplicate native tRNA coordinates')
        categories[key]=cols[note_index]
    return categories


def telomere_calls(aid, records, lengths, table, policy):
    density=float(policy['min_density']);distance=int(policy['max_end_distance']);minimum=int(policy['raw_min_count'])
    if not math.isfinite(density) or not 0<=density<=1 or distance<0:raise ValueError('Invalid telomere density/end distance')
    windows=[];last={};seen=set();window_size=int(policy['window'])
    with Path(table).open() as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            seq=r['id'];end=int(r['window']);start=last.get(seq,0)+1
            if seq not in lengths or end!=min(start-1+window_size,lengths[seq]) or end<start:raise ValueError('Invalid/incomplete tidk window coordinates')
            last[seq]=end;seen.add(seq)
            fwd,rev=int(r['forward_repeat_number']),int(r['reverse_repeat_number']);n=fwd+rev
            if min(fwd,rev)<0 or r['telomeric_repeat']!=policy['motif']:raise ValueError('Invalid tidk counts/motif')
            coverage=n*len(policy['motif'])/(end-start+1)
            if coverage>1+1e-9:raise ValueError('tidk motif coverage exceeds window length')
            if n>=minimum:windows.append(dict(seq=seq,start=start,end=end,fwd=fwd,rev=rev,density=coverage,reason='low_density' if coverage<density else ''))
    if seen!=set(lengths) or any(last[s]!=lengths[s] for s in seen):raise ValueError('Incomplete tidk windows')
    groups=[]
    for i,w in enumerate(windows):
        if w['reason']:continue
        if groups and windows[groups[-1][-1]]['seq']==w['seq'] and windows[groups[-1][-1]]['end']+1==w['start']:groups[-1].append(i)
        else:groups.append([i])
    kept=[];audit=[];rejected=[];outputs={}
    for group in groups:
        first,last_window=windows[group[0]],windows[group[-1]];seq=first['seq'];start,end=first['start'],last_window['end']
        gap=min(start-1,lengths[seq]-end)
        if gap>distance:
            for i in group:windows[i]['reason']='far_from_contig_end'
            continue
        fid=f'{aid}.telomere.filtered{len(kept)+1}'
        fwd=sum(windows[i]['fwd'] for i in group);rev=sum(windows[i]['rev'] for i in group)
        position='spanning' if start-1<=distance and lengths[seq]-end<=distance else '5_prime' if start-1<=distance else '3_prime'
        attrs=dict(ID=fid,Name='telomeric_repeat',motif=policy['motif'],n_repeats=str(fwd+rev),forward_count=str(fwd),reverse_count=str(rev),position=position,
                   motif_density=f'{(fwd+rev)*len(policy["motif"])/(end-start+1):.6f}',end_distance=str(gap))
        kept.append(([seq,'tidk','telomeric_repeat',str(start),str(end),str(fwd+rev),'+' if fwd>=rev else '-','.'],attrs))
        for i in group:outputs[i]=fid
    for i,w in enumerate(windows):
        fid=f'{aid}.telomere.window{i+1}';reason=w['reason'];gap=min(w['start']-1,lengths[w['seq']]-w['end'])
        audit.append(dict(feature_id=fid,decision='rejected' if reason else 'retained',reason=reason or 'passed',output_id=outputs.get(i,''),
                          contig=w['seq'],start=w['start'],end=w['end'],density=f'{w["density"]:.6f}',end_distance=gap,native_category=''))
        if reason:
            rejected.append(([w['seq'],'tidk','telomeric_repeat',str(w['start']),str(w['end']),str(w['fwd']+w['rev']),'.','.'],
                             dict(ID=fid,motif=policy['motif'],qc_rejection=reason)))
    return kept,rejected,audit,'candidate_windows'


def filter_track(aid,track,source,identity,native_status,policy,aux,outdir):
    data,lengths=dictionary(identity);raw=read_gff(source,lengths)
    native=json.loads(Path(native_status).read_text())
    if native.get('status')!='ready' or native['assembly_id']!=aid or native['track']!=track or native['reference_identity']!=data['reference_identity'] or native['gff3_sha256']!=sha(source):raise ValueError('Raw track provenance mismatch')
    audit=[];kept=[];rejected=[];unit='features'
    if not policy['enabled']:
        kept=raw
        audit=[dict(feature_id=a['ID'],decision='retained',reason='filter_disabled',output_id=a['ID'],contig=c[0],start=c[3],end=c[4],density='',end_distance='',native_category='') for c,a in raw]
    elif track=='telomere':
        if len(aux)!=1:raise ValueError('Telomere filtering requires native tidk windows')
        kept,rejected,audit,unit=telomere_calls(aid,raw,lengths,aux[0],policy)
    else:
        categories=trna_categories(aux[0]) if track=='trnascan' and len(aux)==1 else {}
        if track=='trnascan' and len(aux)!=1:raise ValueError('tRNA filtering requires native confidence classifications')
        decisions={};native_tags={}
        for c,a in raw:
            if 'Parent' in a:continue
            reason=''
            if track=='trnascan':
                key=(c[0],int(c[3]),int(c[4]),c[6]);tag=categories.get(key,'not_in_native_confidence_output')
                if key not in categories and c[0] not in ('chrM','M','chrMT','MT'):raise ValueError('Native confidence output lacks tRNA coordinates: '+str(key))
                native_tags[a['ID']]=tag
                if 'high confidence set' not in [v.strip() for v in tag.split(',')]:reason='native_tRNA_not_high_confidence'
            elif track=='infernal':
                if a.get('inc') not in ('!','?'):raise ValueError('Rfam inclusion annotation missing')
                if a.get('olp') not in ('*','~','='):raise ValueError('Rfam overlap annotation missing')
                if a['inc']!='!':reason='below_native_inclusion'
                elif policy['remove_overlaps'] and a['olp']=='=':reason='native_lower_ranked_overlap'
            elif track=='centromere':
                value=float(a['identity'])
                if not math.isfinite(value) or not 0<=value<=1:raise ValueError('Invalid native repeat identity')
                if value<float(policy['min_identity']):reason='low_native_identity'
                elif c[2]=='repeat_region' and policy['high_confidence_hor']:
                    confidence=a.get('confidence','').lower()
                    if confidence not in ('high','low'):raise ValueError('Native HOR confidence missing/unknown')
                    if confidence!='high':reason='native_low_confidence_HOR'
            decisions[a['ID']]=reason
        def reason_for(a,trail=()):
            fid=a['ID']
            if fid in decisions:return decisions[fid]
            if fid in trail:raise ValueError('Cyclic track parents')
            parents=a.get('Parent','').split(',')
            if not all(p in by_id for p in parents):raise ValueError('Unknown track parent')
            reasons=[reason_for(by_id[p],trail+(fid,)) for p in parents]
            if any(reasons) and not all(reasons):raise ValueError('Mixed retained/rejected shared parents require explicit projection')
            return 'parent_rejected' if any(reasons) else ''
        by_id={a['ID']:a for c,a in raw}
        for c,a in raw:
            reason=reason_for(a);(rejected if reason else kept).append((c,a))
            audit.append(dict(feature_id=a['ID'],decision='rejected' if reason else 'retained',reason=reason or 'passed',output_id='' if reason else a['ID'],contig=c[0],start=c[3],end=c[4],density='',end_distance='',native_category=native_tags.get(a['ID'],'')))
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True);details=out/f'{track}_filter';details.mkdir(exist_ok=True)
    selected=out/f'{aid}.{track}.gff3';write_gff(selected,kept);write_gff(details/'rejected.gff3',rejected)
    fields=['feature_id','decision','reason','output_id','contig','start','end','density','end_distance','native_category']
    with (details/'decisions.tsv').open('w',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(audit)
    state=dict(status='ready',assembly_id=aid,track=track,reference_identity=data['reference_identity'],
               features=len(kept),feature_types=dict(Counter(c[2] for c,a in kept)),raw_features=len(raw),raw_feature_types=native['feature_types'],
               source_sha256=sha(source),gff3_sha256=sha(selected),native_source_sha256=native['source_sha256'],filter_policy=policy,
               decision_unit=unit,decision_counts=dict(Counter(r['decision'] for r in audit)),reasons=dict(Counter(r['reason'] for r in audit)),
               auxiliary_sha256=[sha(p) for p in aux])
    dump(out/f'{track}.status.json',state);dump(details/'policy.json',policy)
    return state


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('aid','track','source','identity','native-status','policy','outdir'):p.add_argument('--'+name,required=True)
    p.add_argument('--aux',nargs='*',default=[])
    args=vars(p.parse_args());args['policy']=json.loads(base64.b64decode(args['policy']))
    filter_track(**args)
