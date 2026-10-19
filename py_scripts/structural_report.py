#!/usr/bin/env python3
"""Structural QC adapted from STRUCTURAL_QC and the existing report metrics/plots.

Functional metrics remain unavailable. Structure, AED, junction and repeat flags
remain distinct; no automatic high-confidence list is inferred from them.
"""
import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path
import re
import shutil
import statistics

from collect_annotation_metrics import parse_busco,parse_omark
from model_selection import table,write_table
from prepare_protein_evidence import checksum


def numbers(rows,key):
    values=[float(r[key]) for r in rows if r.get(key) not in (None,'','.')] 
    if any(not math.isfinite(v) for v in values):raise ValueError('Nonfinite structural QC metric: '+key)
    return values


def report(aid,identified,selection,before_table,before_status,before_qc,after_qc,outdir):
    identified,selection,out=Path(identified),Path(selection),Path(outdir)
    ident=json.loads((identified/'annotation_status.json').read_text());sel=json.loads((selection/'selection_status.json').read_text())
    ast=json.loads(Path(before_status).read_text())
    if (ident.get('stage')!='identified_annotation' or ident.get('status')!='ready' or ident['assembly_id']!=aid or
        ast.get('status')!='ready' or ast.get('stage')!='model_assessment' or ast.get('assembly_id')!=aid or
        ident['selection_sha256']!=checksum(selection/'selection_status.json') or sel['input_model_sha256']!=ast['model_sha256'] or
        ast['reference_identity']!=ident['reference_identity'] or checksum(before_table)!=ast['outputs']['transcript_assessment.tsv']):
        raise ValueError('Structural report source identities differ')
    for name,sha in ident['outputs'].items():
        if checksum(identified/name)!=sha:raise ValueError('Identified output changed: '+name)
    for name,sha in sel['outputs'].items():
        if checksum(selection/name)!=sha:raise ValueError('Selected output changed: '+name)
    _,before=table(before_table,'transcript_id');_,after=table(identified/'transcript_assessment.tsv','transcript_id')
    _,decisions=table(selection/'transcript_decisions.tsv','transcript_id')
    if {r['transcript_id'] for r in before}!={r['transcript_id'] for r in decisions}:raise ValueError('Decision/assessment IDs differ')
    expected={r['transcript_id'] for r in decisions if r['decision']=='retained'}
    if expected!={r['original_transcript_id'] for r in after}:raise ValueError('Final assessment differs from selection')
    out.mkdir(parents=True,exist_ok=True);native=out/'native';native.mkdir(exist_ok=True)
    qc_metrics=[];qc_inputs={};raw=[]
    for phase,paths,expected_sha in [('before',before_qc,sel['input_model_sha256']),('after',after_qc,ident['gff3_sha256'])]:
        files={Path(p).name:Path(p) for p in paths};qc_inputs[phase]={name:checksum(p) for name,p in files.items()}
        metadata=files.get('selection.json')
        if metadata:
            meta=json.loads(metadata.read_text())
            if meta['assembly_id']!=aid or meta['gff3_sha256']!=expected_sha:raise ValueError('QC assessed a different annotation')
            if phase=='after' and meta['coding_transcripts']!=len(after):raise ValueError('QC transcript count differs')
            for tool,suffix,parser in [('busco',f'{aid}_busco.txt',parse_busco),('omark',f'{aid}.summary.txt',parse_omark)]:
                path=files.get(suffix);enabled=meta['settings'].get('run_'+tool,False)
                if enabled and not path:raise ValueError('Missing enabled '+phase+' '+tool+' QC result')
                if path:
                    text=path.read_text()
                    if tool=='busco' and not re.search(r'C:[\d.]+%.*S:[\d.]+%.*D:[\d.]+%.*F:[\d.]+%.*M:[\d.]+%.*n:\s*\d+',text):raise ValueError('Unrecognized BUSCO summary')
                    if tool=='omark' and (not re.search(r'^S:[\d.]+%.*D:.*M:',text,re.M) or not re.search(r'^A:[\d.]+%.*I:.*C:.*U:',text,re.M)):raise ValueError('Unrecognized OMArk summary')
                    values=parser(path)
                    if tool=='omark':values={k:v for k,v in values.items() if k in ('omark_single','omark_duplicated','omark_missing','omark_consistent','omark_inconsistent','omark_contaminant','omark_unknown')}
                    for key,value in values.items():qc_metrics.append(dict(model_set=phase,metric=key,value=value))
            if f'{aid}_annotation_stats' not in files:raise ValueError('Missing '+phase+' AGAT statistics')
        elif files or phase=='after':raise ValueError('QC files require representative-selection provenance')
        for name,path in files.items():
            dest=native/f'{phase}_{name}';shutil.copyfile(path,dest)
            if name.endswith(('_busco.txt','.summary.txt','_detailed.txt','_annotation_stats')):
                raw.append((phase+' '+name,dest.name,path.read_text()))
    # Never replace unavailable native tool metrics by zeros.
    write_table(out/'completeness_metrics.tsv',['model_set','metric','value'],qc_metrics)
    assembly_rows=[];gene_rows=[]
    for label,rows in [('before',before),('after',after)]:
        groups=defaultdict(list)
        for r in rows:groups[r['gene_id']].append(r)
        summary=dict(model_set=label,genes=len(groups),transcripts=len(rows),single_exon_transcripts=sum(r['is_single_exon']=='Y' for r in rows),
            structure_pass=sum(r['structure_flag']=='PASS' for r in rows),structure_warn=sum(r['structure_flag']=='WARN' for r in rows),structure_fail=sum(r['structure_flag']=='FAIL' for r in rows),
            te_flagged_genes=sum(any(r['te_gene_flag']=='True' for r in rs) for rs in groups.values()) if any(r['te_gene_flag']!='.' for r in rows) else '.',
            functional_support='pending',protein_length_median=statistics.median(numbers(rows,'protein_length_aa')),
            exons_per_transcript_median=statistics.median(numbers(rows,'n_exons')))
        assembly_rows.append(summary)
        if label=='after':
            rank={'PASS':0,'WARN':1,'FAIL':2}
            for gid,rs in sorted(groups.items()):
                flags=[r['structure_flag'] for r in rs]
                row=dict(gene_id=gid,original_gene_id=rs[0]['original_gene_id'],transcripts=len(rs),
                    best_structure_flag=min(flags,key=rank.get),worst_structure_flag=max(flags,key=rank.get),
                    te_flagged=any(r['te_gene_flag']=='True' for r in rs) if any(r['te_gene_flag']!='.' for r in rs) else '.',functional_support='pending')
                for metric in ('rna_aed','protein_aed','junction_coding_pass_fraction'):
                    v=numbers(rs,metric);row[metric+'_min']=min(v) if v else '.';row[metric+'_max']=max(v) if v else '.'
                gene_rows.append(row)
    write_table(out/'structural_qc.tsv',list(assembly_rows[0]),assembly_rows)
    write_table(out/'structural_qc_genes.tsv',list(gene_rows[0]),gene_rows)
    shutil.copyfile(identified/'transcript_assessment.tsv',out/'structural_qc_transcripts.tsv')
    # Reuse the original report's before/after distributions, consistent colors,
    # separate completeness panels and gene best/worst aggregation.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, NullLocator, FuncFormatter
    import numpy as np
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    plots=out/'plots';plots.mkdir(exist_ok=True);figures=[];colors=['#466A95','#27856C']
    def save(fig,name,title,caption):
        fig.suptitle(title);fig.tight_layout(rect=(0,0,1,.95))
        for ext in ('png','svg'):fig.savefig(plots/f'{name}.{ext}',dpi=130,bbox_inches='tight')
        plt.close(fig);figures.append((name,title,caption))
    def ecdf(ax,v,label,color):
        v=np.sort(v)
        if len(v):ax.step(np.r_[v[0],v],np.r_[0,np.arange(1,len(v)+1)/len(v)],where='post',label=label,color=color)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,key,title in zip(axes,['genes','transcripts'],['Genes','Coding transcripts']):
        values=[r[key] for r in assembly_rows];ax.bar(['Before','After'],values,color=colors);ax.set_title(title);ax.set_ylabel('Count')
        for i,v in enumerate(values):ax.text(i,v,f'{v:,}',ha='center',va='bottom')
        ax.margins(y=.2);ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    save(fig,'retention','Model counts before and after selection',f"{sel['excluded_transcripts']:,} transcripts and {sel['excluded_genes']:,} genes excluded; a gene survives if any isoform survives.")
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,key,title in zip(axes,['protein_length_aa','cds_length_bp','n_exons'],['Protein length (aa)','CDS length per transcript (bp)','Exons per transcript']):
        for rows,label,color in zip([before,after],['Before','After'],colors):
            ecdf(ax,numbers(rows,key),label,color)
        ax.set_xlabel(title);ax.set_ylabel('Fraction of transcripts ≤ value');ax.set_ylim(0,1);ax.legend()
        if key=='n_exons':ax.xaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
        else:
            v=numbers(before,key);ax.set_xscale('log');ax.set_xticks(sorted(set(round(x) for x in np.geomspace(min(v),max(v),4))))
            ax.xaxis.set_minor_locator(NullLocator());ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'))
    save(fig,'structure','Structure distributions before and after selection','Cumulative distributions use all coding transcripts; logarithmic length axes keep long tails visible. Exon counts use a linear axis.')
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,key,title in zip(axes,['n_introns','min_intron_len','isoforms'],['Introns per transcript','Shortest intron per spliced transcript (bp)','Transcripts per gene']):
        for rows,label,color in zip([before,after],['Before','After'],colors):
            v=list(Counter(r['gene_id'] for r in rows).values()) if key=='isoforms' else numbers(rows,key)
            ecdf(ax,v,label,color)
        ax.set(xlabel=title,ylabel='Fraction ≤ value',ylim=(0,1));ax.xaxis.set_major_locator(MaxNLocator(nbins=4,integer=True))
        if ax.lines:ax.legend()
        else:ax.text(.5,.5,'No applicable introns',ha='center',transform=ax.transAxes)
    save(fig,'complexity','Intron and isoform structure','Intron counts include single-exon transcripts at zero. Minimum intron length excludes single-exon transcripts. Isoform counts use genes as their denominator.')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,key,title in zip(axes,['rna_aed','protein_aed'],['RNA AED','Protein AED']):
        for rows,label,color in zip([before,after],['Before','After'],colors):
            v=np.sort(numbers(rows,key))
            if len(v):ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=f'{label} (n={len(v):,})',color=color)
        ax.set(xlabel=title+' (lower is closer agreement)',ylabel='Fraction of assessed transcripts ≤ AED',xlim=(0,1),ylim=(0,1))
        if ax.lines:ax.legend()
        else:ax.text(.5,.5,'Evidence unavailable',ha='center',transform=ax.transAxes)
    save(fig,'aed','Separate RNA and protein evidence agreement','Missing evidence is excluded from the assessed denominator; it is not assigned AED=1. RNA and protein scores are not combined.')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,kind in zip(axes,['coding','utr']):
        for rows,label,color in zip([before,after],['Before','After'],colors):
            v=np.sort(numbers(rows,'junction_'+kind+'_pass_fraction'))
            if len(v):ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=f'{label} (n={len(v):,})',color=color)
        ax.set(xlabel=f'{kind.upper()} junction pass fraction',ylabel='Fraction of applicable transcripts ≤ value',xlim=(0,1),ylim=(0,1))
        if ax.lines:ax.legend()
        else:ax.text(.5,.5,'Not applicable / no RNA',ha='center',transform=ax.transAxes)
    save(fig,'junctions','Coding and UTR-associated splice support','Only assessed, applicable transcript fractions are plotted. A shared genomic junction may occur in several transcripts; these are transcript summaries, not counts of unique genomic junctions.')
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,kind in zip(axes,['coding','utr']):
        cats=['pass','fail','not_observed','unresolved_strand'];x=np.arange(len(cats));available=False
        for i,rows in enumerate([before,after]):
            if not any(r.get('junction_'+kind+'_pass') not in ('.',None,'') for r in rows):continue
            counts=[sum(numbers(rows,'junction_'+kind+'_'+c)) for c in cats]
            ax.bar(x+(i-.5)*.36,counts,.36,label=['Before','After'][i],color=colors[i]);available=True
        ax.set_xticks(x,['Pass','Fail','Not observed','Unresolved strand'],rotation=20,ha='right')
        ax.set_title(kind.upper()+' junction occurrences');ax.set_ylabel('Transcript junction occurrences');ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if available:ax.legend()
        else:ax.text(.5,.5,'RNA support unavailable',ha='center',transform=ax.transAxes)
    save(fig,'junction_occurrences','Splice-junction outcomes','These are occurrences across transcripts, not distinct genomic junctions or read counts. Shared isoforms repeat a junction in this population. Distinct-junction and readset-level results remain in the native Portcullis report.')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,key,title in zip(axes,['rna_aed','protein_aed'],['RNA AED','Protein AED']):
        pairs=[(float(r[key]),float(r['junction_coding_pass_fraction'])) for r in after if r[key]!='.' and r['junction_coding_pass_fraction']!='.']
        if pairs:
            x,y=zip(*pairs);h=ax.hist2d(x,y,bins=30,range=[[0,1],[0,1]],cmap='Blues',cmin=1);fig.colorbar(h[3],ax=ax,label='Transcripts per bin')
        else:ax.text(.5,.5,'Paired evidence unavailable',ha='center',transform=ax.transAxes)
        ax.set(xlabel=title,ylabel='Coding junction pass fraction',xlim=(0,1),ylim=(0,1))
    save(fig,'evidence_agreement','AED and splice support provide different evidence','Selected models only. Density bins expose concentrations without overplotting. Neither axis alone establishes biological correctness.')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,keys,title in [(axes[0],['PASS','WARN','FAIL'],'Transcript structural flags'),(axes[1],['repeat_overlap','repeat_library_homology','repeat_overlap;repeat_library_homology'],'Gene repeat-association flags')]:
        x=np.arange(len(keys))
        for i,rows in enumerate([before,after]):
            if ax is axes[0]:counter=Counter(r['structure_flag'] for r in rows)
            else:
                if not any(r['te_gene_flag']!='.' for r in rows):continue
                counter=Counter({r['gene_id']:r['te_evidence'] for r in rows}.values())
            ax.bar(x+(i-.5)*.36,[counter[k] for k in keys],.36,label=['Before','After'][i],color=colors[i])
        labels=keys if ax is axes[0] else ['Overlap','Homology','Both']
        ax.set_xticks(x,labels);ax.set_ylabel('Transcripts' if ax is axes[0] else 'Genes');ax.set_title(title);ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if ax.containers:ax.legend()
        else:ax.text(.5,.5,'TE assessment unavailable',ha='center',transform=ax.transAxes)
    save(fig,'flags','Diagnostic flags retained for interpretation','PASS/WARN/FAIL are diagnostic severity. Repeat flags are inherited from the assessed unfiltered locus and can reflect an excluded isoform; overlap and library homology are broad association signals, not automatic confidence labels.')
    fig,ax=plt.subplots(figsize=(9,4));reasons=sel['removal_reason_counts']
    reason_labels={'single_exon_short_cds':f"Single-exon CDS < {sel['settings']['min_cds_length']:g} bp",'short_gene_span':f"Gene span < {sel['settings']['min_gene_length']:g} bp",'all_introns_short':f"All introns < {sel['settings']['min_intron_length']:g} bp",'gene_repeat_association':'Gene repeat association','high_rna_aed':'RNA AED above threshold','high_protein_aed':'Protein AED above threshold','low_coding_junction_support':'Coding junction support below threshold'}
    if reasons:ax.barh([reason_labels.get(k,k) for k in reasons],list(reasons.values()),color=colors[0]);ax.set_xlabel('Excluded transcripts triggering the rule')
    else:ax.text(.5,.5,'No transcripts excluded',ha='center',transform=ax.transAxes)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    save(fig,'decisions','Explicit selection reasons','A transcript can trigger more than one rule; reason counts need not sum to the total excluded. Full decisions and final ID maps retain all input models.')
    fig,axes=plt.subplots(1,3,figsize=(16,4))
    for ax,tool,keys in [(axes[0],'BUSCO',['busco_single','busco_duplicated','busco_fragmented','busco_missing']),
                          (axes[1],'OMArk completeness',['omark_single','omark_duplicated','omark_missing']),
                          (axes[2],'OMArk consistency',['omark_consistent','omark_inconsistent','omark_contaminant','omark_unknown'])]:
        x=np.arange(len(keys));available=False
        for i,phase in enumerate(['before','after']):
            d={r['metric']:r['value'] for r in qc_metrics if r['model_set']==phase}
            if all(k in d for k in keys):ax.bar(x+(i-.5)*.36,[d[k] for k in keys],.36,label=phase.title(),color=colors[i]);available=True
        ax.set_xticks(x,[k.split('_',1)[1].title() for k in keys],rotation=25,ha='right');ax.set_ylabel('%');ax.set_ylim(0,100);ax.set_title(tool)
        if available:ax.legend()
        else:ax.text(.5,.5,'Not requested / unavailable',ha='center',transform=ax.transAxes)
    save(fig,'completeness','Completeness and proteome consistency','BUSCO/OMArk assess one longest protein per coding gene, not every isoform. Percentages use each native tool’s denominator; native summaries and representative tables are supplied.')
    # Native OMArk may omit very short proteins: preserve unclassified explicitly.
    classes=defaultdict(list);ump=next((Path(p) for p in after_qc if Path(p).name==f'{aid}.ump'),None)
    if ump:
        category=None
        for line in ump.read_text().splitlines():
            if line.startswith('>'):category=line[1:].strip()
            elif line.strip() and not line.startswith('#'):
                if not category:raise ValueError('OMArk classes lack a category')
                classes[line.strip()].append(category)
        if not set(classes)<={r['transcript_id'] for r in after}:raise ValueError('OMArk classified unknown final IDs')
    fields=list(after[0])+['omark_representative_classification']
    repfile=next((Path(p) for p in after_qc if Path(p).name=='representatives.tsv'),None)
    reps={r['transcript_id'] for r in table(repfile,'transcript_id')[1]} if repfile else set()
    if not set(classes)<=reps:raise ValueError('OMArk classes differ from assessed representative IDs')
    for r in after:r['omark_representative_classification']=';'.join(classes[r['transcript_id']]) if r['transcript_id'] in classes else ('unclassified' if ump and r['transcript_id'] in reps else 'not_assessed')
    write_table(out/'structural_qc_transcripts.tsv',fields,after)
    summary=dict(status='ready',stage='structural_qc',assembly_id=aid,reference_identity=ident['reference_identity'],
        identified_model_sha256=ident['gff3_sha256'],selection_sha256=ident['selection_sha256'],model_sets=assembly_rows,
        native_qc_sources=qc_inputs,functional_annotation='pending',confidence_list_created=False,
        scope='Selection and final-ID integrity, structure, separate native AED, junction/repeat diagnostics and before/after AGAT/BUSCO/OMArk; no further model removal',
        outputs={str(p.relative_to(out)).replace('\\','/'):checksum(p) for p in out.rglob('*') if p.is_file() and p.name not in ('structural_qc_status.json','README.md')})
    (out/'structural_qc_status.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Structural annotation QC','',f'Assembly: **{aid}**. Status: **ready**.','',
        '| Model set | Genes | Coding transcripts | PASS | WARN | FAIL |','|---|---:|---:|---:|---:|---:|']
    for r in assembly_rows:lines.append('| '+ ' | '.join(str(r[k]) for k in ['model_set','genes','transcripts','structure_pass','structure_warn','structure_fail'])+' |')
    intro='This report compares accepted unfiltered post-UTR models with the selected final-ID set. Retained AED and junction results are re-keyed through the ID map, not recalculated. It records selection effects; it does not establish universal biological accuracy. Functional annotation remains pending. No combined AED or automatic high-confidence list is generated.'
    lines+=['',intro,'','[Model summaries](structural_qc.tsv) | [Gene table](structural_qc_genes.tsv) | [Transcript table](structural_qc_transcripts.tsv) | [Completeness metrics](completeness_metrics.tsv)', '',
        f'[Final identified models and maps](../../../structural_annotation/identified/{aid}/README.md) | [Selection decisions](../../../structural_annotation/selection/{aid}/README.md) | [Before QC](../../final/README.md) | [Post-filter QC](../../selected/README.md) | [Native Portcullis QC](../../junctions/{aid}/README.md)',
        f'[Before AGAT distribution figures](../../final/{aid}/agat/annotation_stats_distribution_plots/{aid}/) | [After AGAT distribution figures](../../selected/{aid}/agat/annotation_stats_distribution_plots/{aid}/)','']
    for name,title,caption in figures:
        lines+=['## '+title,'',caption,'',f'![{title}](plots/{name}.png)',f'[Vector figure](plots/{name}.svg)','']
    for title,name,text in raw:
        lines += [f'[{title}](native/{name})']
    lines+=['','Native AGAT reports include all features and isoforms; native distribution figures remain in the AGAT QC directories. The figures above use explicitly labelled transcript or gene populations. Missing evidence and non-applicable single-exon junction metrics remain distinct from zero support.','', '[Machine-readable validation](structural_qc_status.json)']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('assembly-id','identified','selection','before-table','before-status','outdir'):p.add_argument('--'+k,required=True)
    for k in ('before-qc','after-qc'):p.add_argument('--'+k,nargs='*',default=[])
    a=p.parse_args();report(a.assembly_id,a.identified,a.selection,a.before_table,a.before_status,a.before_qc,a.after_qc,a.outdir)
