#!/usr/bin/env python3
"""Select transcripts using the original structural cutoffs and opt-in evidence rules.

Explicit parents replace the old current-gene/name fallback. Alternative isoforms
are evaluated separately; a gene remains if any transcript remains.
"""
import argparse
import base64
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
from urllib.parse import quote, unquote

from evidence_aed import models
from model_metrics import sequences
from prepare_protein_evidence import checksum


DEFAULTS = dict(structural_action='remove', min_cds_length=300, min_gene_length=200,
    min_intron_length=10, aed_action='remove', aed_source=None, max_aed=None,
    junction_action='remove', min_junction_fraction=None, te_action='flag')


def settings(values):
    unknown = set(values)-set(DEFAULTS)
    if unknown: raise ValueError('Unknown selection settings: '+str(sorted(unknown)))
    result = dict(DEFAULTS, **values)
    for key in ('structural_action','aed_action','junction_action','te_action'):
        if result[key] not in ('off','flag','remove'): raise ValueError('Invalid rule action: '+key)
    for key in ('min_cds_length','min_gene_length','min_intron_length','max_aed','min_junction_fraction'):
        v=result[key]
        if v is None:
            if key.startswith('min_') and key!='min_junction_fraction': raise ValueError('Structural cutoffs require a number; use 0 to disable')
            continue
        v=float(v)
        if not math.isfinite(v) or v<0: raise ValueError('Invalid selection cutoff: '+key)
        if key in ('max_aed','min_junction_fraction'):
            if v>1: raise ValueError('Evidence cutoff must be between 0 and 1')
        elif not v.is_integer(): raise ValueError('Length cutoffs must be whole base pairs')
        result[key]=v
    if result['max_aed'] is not None and result['aed_action']!='off' and result['aed_source'] not in ('rna','protein'):
        raise ValueError('Set filter_aed_source to rna or protein when AED filtering is enabled')
    return result


def table(path,key):
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t'); fields=reader.fieldnames
        rows=list(reader)
    if not fields or any(None in r or any(v is None for v in r.values()) for r in rows): raise ValueError('Malformed assessment table')
    if any(not r.get(key) for r in rows) or len({r[key] for r in rows})!=len(rows): raise ValueError('Duplicate/missing assessment IDs')
    return fields,rows


def write_table(path,fields,rows):
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def gff_partition(source,destination,tx,selected):
    """Retain descendants by explicit ancestry; split shared Parents per output.

    Preserve original gene bounds and feature IDs. A gene may appear in both
    partitions when its alternative isoforms have different decisions.
    """
    lines=Path(source).read_text().splitlines(keepends=True); nodes={}; parsed=[]
    for line in lines:
        if line.startswith('##FASTA') or line.startswith('>'): raise ValueError('Embedded FASTA is not supported')
        if line.startswith('#') or not line.strip(): parsed.append((line,None,None)); continue
        cols=line.rstrip('\n\r').split('\t')
        if len(cols)!=9: raise ValueError('Expected nine GFF columns')
        attrs={k.strip():v for k,v in (part.split('=',1) for part in cols[8].split(';') if '=' in part)}
        fid=unquote(attrs['ID']) if 'ID' in attrs else None
        parents=tuple(unquote(p) for p in attrs.get('Parent','').split(',') if p)
        node=(cols[2],parents)
        if fid and fid in nodes and nodes[fid]!=node: raise ValueError('Inconsistent multipart GFF feature ID')
        if fid:nodes[fid]=node
        if cols[2]!='gene' and not parents: raise ValueError('Unparented GFF feature')
        parsed.append((line,cols,(fid,parents)))
    kept_genes={tx[t]['gene_id'] for t in selected}; memo={}; visiting=set()
    def retain(fid):
        if fid in memo:return memo[fid]
        if fid in visiting or fid not in nodes: raise ValueError('Cyclic/unresolved GFF parent')
        visiting.add(fid);kind,parents=nodes[fid]
        if kind=='gene':result=fid in kept_genes
        elif kind in ('mRNA','transcript'):result=fid in selected
        else:result=any([retain(p) for p in parents])
        visiting.remove(fid);memo[fid]=result;return result
    # Resolve every node even if excluded, so malformed ancestry cannot disappear.
    for fid in nodes:retain(fid)
    with Path(destination).open('w') as handle:
        for line,cols,info in parsed:
            if cols is None:handle.write(line);continue
            fid,parents=info
            chosen=tuple(p for p in parents if retain(p))
            keep=retain(fid) if fid else bool(chosen)
            if not keep:continue
            if parents and chosen!=parents:
                attrs=cols[8].split(';')
                cols[8]=';'.join('Parent='+','.join(quote(p,safe='.:_-') for p in chosen) if a.split('=',1)[0].strip()=='Parent' else a for a in attrs)
                handle.write('\t'.join(cols)+'\n')
            else:handle.write(line)


def select(aid,gff,proteins,cds,validation,identity,assessment,assessment_status,outdir,options=None):
    opts=settings(options or {})
    state=json.loads(Path(validation).read_text()); ast=json.loads(Path(assessment_status).read_text())
    ref=json.loads(Path(identity).read_text()); lengths={r['id']:r['length'] for r in ref['contigs']}
    if (state.get('stage')!='final_annotation' or ast.get('stage')!='model_assessment' or
        any(s.get('status')!='ready' or s.get('assembly_id')!=aid or s.get('reference_identity')!=ref['reference_identity'] for s in (state,ast))):
        raise ValueError('Selection requires accepted annotation/assessment/reference')
    if checksum(gff)!=state['gff3_sha256'] or ast['model_sha256']!=state['gff3_sha256']: raise ValueError('Assessment and model identities differ')
    if checksum(assessment)!=ast['outputs']['transcript_assessment.tsv']: raise ValueError('Assessment table changed')
    for path,key in [(proteins,'proteins_sha256'),(cds,'cds_sha256')]:
        if checksum(path)!=state[key]: raise ValueError('Accepted coding sequences changed')
    fields,rows=table(assessment,'transcript_id'); tx=models(gff,lengths)
    if set(tx)!={r['transcript_id'] for r in rows}: raise ValueError('Assessment IDs differ from models')
    protein,coding=sequences(proteins),sequences(cds)
    if set(tx)!=set(protein) or set(tx)!=set(coding): raise ValueError('Model/protein/CDS IDs differ')
    aed_on=opts['max_aed'] is not None and opts['aed_action']!='off'
    junction_on=opts['min_junction_fraction'] is not None and opts['junction_action']!='off'
    if aed_on and not ast['aed_enabled']:raise ValueError('AED selection requires enabled AED assessment')
    if junction_on and not ast['junction_qc_enabled']:raise ValueError('Junction selection requires enabled junction QC')
    if opts['te_action']=='remove' and not ast['te_enabled']:raise ValueError('TE removal requires enabled TE assessment')
    decisions=[];kept=set();by_gene=defaultdict(list)
    for row in rows:
        tid=row['transcript_id']; m=tx[tid]
        if row['gene_id']!=m['gene_id']: raise ValueError('Assessment gene parent differs')
        exons=sorted(set(m['exons'])); introns=[r[0]-l[1] for l,r in zip(exons,exons[1:])]
        cds_bp=sum(e-s for s,e in set(m['cds']))
        if int(row['cds_length_bp'])!=cds_bp or int(row['n_introns'])!=len(introns):raise ValueError('Assessment structure differs')
        reasons=[];remove=[];skipped=[]
        def hit(reason,action):
            reasons.append(reason+':'+action)
            if action=='remove':remove.append(reason)
        if opts['structural_action']!='off':
            if len(exons)==1 and 0<cds_bp<opts['min_cds_length']:hit('single_exon_short_cds',opts['structural_action'])
            if int(row['gene_span_bp'])<opts['min_gene_length']:hit('short_gene_span',opts['structural_action'])
            if introns and max(introns)<opts['min_intron_length']:hit('all_introns_short',opts['structural_action'])
        if aed_on:
            value=row[opts['aed_source']+'_aed']
            if value=='.':skipped.append(opts['aed_source']+'_aed_unavailable')
            else:
                value=float(value)
                if not math.isfinite(value) or not 0<=value<=1:raise ValueError('Invalid AED score')
                if value>opts['max_aed']:hit('high_'+opts['aed_source']+'_aed',opts['aed_action'])
        if junction_on:
            value=row['junction_coding_pass_fraction']
            if value=='.':skipped.append('coding_junctions_not_applicable' if int(row['junction_coding_junctions'])==0 else 'coding_junction_support_unavailable')
            else:
                value=float(value)
                if not math.isfinite(value) or not 0<=value<=1:raise ValueError('Invalid coding junction fraction')
                if value<opts['min_junction_fraction']:hit('low_coding_junction_support',opts['junction_action'])
        if opts['te_action']!='off':
            if row['te_gene_flag']=='.':skipped.append('te_assessment_unavailable')
            elif row['te_gene_flag']=='True':hit('gene_repeat_association',opts['te_action'])
            elif row['te_gene_flag']!='False':raise ValueError('Invalid gene TE flag')
        decision=dict(gene_id=m['gene_id'],transcript_id=tid,decision='excluded' if remove else 'retained',
            removal_reasons=';'.join(remove) or '.',rule_hits=';'.join(reasons) or '.',
            unavailable_or_not_applicable=';'.join(skipped) or '.',structure_flag=row['structure_flag'])
        decisions.append(decision);by_gene[m['gene_id']].append(decision)
        if not remove:kept.add(tid)
    excluded=set(tx)-kept; out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    for label,ids in [('retained',kept),('excluded',excluded)]:
        folder=out/label;folder.mkdir(exist_ok=True)
        gff_partition(gff,folder/f'{aid}.gff3',tx,ids)
        parsed=models(folder/f'{aid}.gff3',lengths) if ids else {}
        if set(parsed)!=ids or any(parsed[t]['cds']!=tx[t]['cds'] or parsed[t]['exons']!=tx[t]['exons'] for t in ids):raise ValueError('Selection altered retained structure')
        for name,seqs in [('proteins',protein),('cds',coding)]:
            with (folder/f'{aid}_{name}.fa').open('w') as handle:
                for tid,seq in seqs.items():
                    if tid in ids:handle.write(f'>{tid}\n{seq}\n')
        write_table(folder/'transcript_assessment.tsv',fields,[r for r in rows if r['transcript_id'] in ids])
    write_table(out/'transcript_decisions.tsv',list(decisions[0]),decisions)
    gene_rows=[dict(gene_id=gid,transcripts=len(ds),retained_transcripts=sum(d['decision']=='retained' for d in ds),
        excluded_transcripts=sum(d['decision']=='excluded' for d in ds),
        decision='retained' if any(d['decision']=='retained' for d in ds) else 'excluded',
        removal_reasons=';'.join(sorted({reason for d in ds for reason in d['removal_reasons'].split(';') if reason!='.'})) or '.') for gid,ds in sorted(by_gene.items())]
    write_table(out/'gene_decisions.tsv',list(gene_rows[0]),gene_rows)
    result=dict(status='ready',stage='model_selection',assembly_id=aid,reference_identity=ref['reference_identity'],
        input_model_sha256=checksum(gff),assessment_status_sha256=checksum(assessment_status),settings=opts,
        input_transcripts=len(tx),retained_transcripts=len(kept),excluded_transcripts=len(excluded),
        input_genes=len(by_gene),retained_genes=sum(r['decision']=='retained' for r in gene_rows),
        excluded_genes=sum(r['decision']=='excluded' for r in gene_rows),empty_retained_set=not bool(kept),
        removal_reason_counts=dict(Counter(v for r in decisions for v in r['removal_reasons'].split(';') if v!='.')),
        scope='Original structural cutoffs per transcript; gene retained if any isoform survives; no renaming or coding changes',
        outputs={str(p.relative_to(out)).replace('\\','/'):checksum(p) for p in [
            *[out/label/name for label in ('retained','excluded') for name in (f'{aid}.gff3',f'{aid}_proteins.fa',f'{aid}_cds.fa','transcript_assessment.tsv')],
            out/'transcript_decisions.tsv',out/'gene_decisions.tsv']})
    (out/'selection_status.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'README.md').write_text('# Model selection\n\nStatus: **ready**\n\n'
        f"Transcripts: {len(tx)} input, **{len(kept)} retained**, {len(excluded)} excluded. "
        f"Genes: {len(by_gene)} input, **{result['retained_genes']} retained**, {result['excluded_genes']} excluded.\n\n"
        '[Transcript decisions](transcript_decisions.tsv) | [Gene decisions](gene_decisions.tsv) | [Settings and validation](selection_status.json)\n\n'
        f'[Retained GFF](retained/{aid}.gff3) | [Proteins](retained/{aid}_proteins.fa) | [CDS](retained/{aid}_cds.fa) | '
        '[Retained assessment](retained/transcript_assessment.tsv)\n\n'
        f'[Excluded GFF](excluded/{aid}.gff3) | [Proteins](excluded/{aid}_proteins.fa) | [CDS](excluded/{aid}_cds.fa)\n\n'
        'Selection evaluates each isoform before deciding whether to exclude a gene. '
        'Shared feature parents are split between outputs; original IDs, gene bounds, coordinates and coding sequences are preserved. '
        'All triggered rules are reported. Diagnostic structure FAIL/WARN is not itself a removal rule. '
        'Missing evidence is not a score of 1 or zero junction support: the rule is explicitly skipped. '
        'Junction selection uses coding junctions only; single-exon/no-coding-junction models are not excluded by that rule. '
        'TE removal, if explicitly enabled, applies a gene-level flag to all its isoforms.\n\n'
        'Original accepted annotations and evidence reports remain preserved. Empty partitions are valid; an empty retained set is identified in the status. '
        'Final identifiers, selected-set QC and functional annotation follow as separate stages.\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','gff','proteins','cds','validation','identity','assessment','assessment-status','outdir','settings'):
        p.add_argument('--'+key,required=True)
    a=p.parse_args()
    select(a.assembly_id,a.gff,a.proteins,a.cds,a.validation,a.identity,a.assessment,a.assessment_status,a.outdir,
        json.loads(base64.b64decode(a.settings,validate=True)))
