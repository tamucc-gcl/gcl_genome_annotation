#!/usr/bin/env python3
"""Results/methods-only Markdown report from accepted native stage outputs."""
import argparse
import base64
from collections import Counter
import csv
import json
import math
import re
from pathlib import Path
import statistics
import textwrap
from urllib.parse import unquote

from generate_annotation_report import fmt, table
from reference_repeats import sha, dump


def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def escape(value): return str(value).replace('|', '\\|').replace('\r', ' ').replace('\n', ' ')
def display(value):
    if value is None or value == '.': return '—'
    if isinstance(value, bool): return 'yes' if value else 'no'
    if isinstance(value, int): return fmt(value)
    if isinstance(value, float):
        if not math.isfinite(value): raise ValueError('Nonfinite report metric')
        return f'{value:.6g}'
    return escape(value)
def md_table(headers, rows): return table(headers, [[display(v) for v in r] for r in rows])
def tsv(path):
    with Path(path).open(encoding='utf-8') as h: return list(csv.DictReader(h, delimiter='\t'))
def by_assembly(state, aid):
    rows = [r for r in state.get('assemblies', []) if r['assembly_id'] == aid]
    if len(rows) != 1: raise ValueError('Missing/duplicate assembly in stage report: '+aid)
    return rows[0]


def build(request, outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    assets = out/'annotation_report_plots'; assets.mkdir(exist_ok=True)
    data = out/'annotation_report_data'; data.mkdir(exist_ok=True)
    stages = {}; provenance = {}
    for key, aid, path in request['stages']:
        if (key, aid) in stages: raise ValueError('Duplicate report stage: '+key+'/'+aid)
        state = read(path)
        if state.get('status') not in ('ready', 'not_requested'): raise ValueError('Unaccepted report stage: '+key+'/'+aid)
        stages[key, aid] = state; provenance[key+'/'+aid] = sha(path)
        dump(data/f'{key}.{aid}.json', state)
    for key in ('preflight','protein','rna','rna_alignment','reference_repeats','braker'):
        if (key, 'run') not in stages: raise ValueError('Missing report stage: '+key)
    preflight = stages['preflight','run']
    expected = sorted(r['assembly_id'] for r in preflight['assemblies'])
    releases = {Path(p).name: Path(p) for p in request['releases']}
    if len(releases) != len(request['releases']) or sorted(releases) != expected: raise ValueError('Report releases differ from requested assemblies')
    structural = {aid: Path(p) for aid, p in request['structural']}
    if len(structural) != len(request['structural']) or sorted(structural) != expected: raise ValueError('Missing/duplicate structural report')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'svg.fonttype':'none'})
    figures = []; metrics = []; parts = []; overview = []
    def save(fig, aid, name):
        stem = aid+'.'+name; fig.tight_layout()
        for ext in ('png','svg'): fig.savefig(assets/(stem+'.'+ext), dpi=160, bbox_inches='tight')
        plt.close(fig); figures.append(stem)
        return f'![{escape(aid)}: {name}](reports/annotation_report_plots/{stem}.png)\n\n[Vector figure](reports/annotation_report_plots/{stem}.svg)'
    for ordinal, aid in enumerate(expected, 1):
        folder = releases[aid]; release = read(folder/'release_status.json')
        if release.get('stage') != 'annotation_release' or release.get('status') != 'ready' or release['assembly_id'] != aid:
            raise ValueError('Report requires accepted release: '+aid)
        ref = release['reference_identity']; provenance['release/'+aid] = sha(folder/'release_status.json')
        def payload(name):
            p = folder/name
            if name not in release['outputs'] or sha(p) != release['outputs'][name]['sha256']: raise ValueError('Release report input changed: '+name)
            return p
        ident = read(payload('provenance/annotation_status.json'))
        track = read(payload('provenance/genome_tracks_status.json'))
        repeat = read(payload('provenance/reference_repeat_status.json'))
        reference = read(payload('provenance/reference.json'))
        if any(s['reference_identity'] != ref for s in (ident, track, repeat, reference)): raise ValueError('Report reference mismatch')
        sq = structural[aid]; sqs = read(sq/'structural_qc_status.json')
        if sqs.get('status') != 'ready' or sqs['assembly_id'] != aid or sqs['reference_identity'] != ref or sqs['identified_model_sha256'] != ident['gff3_sha256']:
            raise ValueError('Structural report assessed a different release')
        selected_sets=[r for r in sqs['model_sets'] if r['model_set']=='after']
        if len(selected_sets)!=1 or (selected_sets[0]['genes'],selected_sets[0]['transcripts'])!=(release['genes'],release['coding_transcripts']):
            raise ValueError('Structural report counts differ from release')
        def qcfile(name):
            p = sq/name
            if name not in sqs['outputs'] or sha(p) != sqs['outputs'][name]: raise ValueError('Structural report input changed: '+name)
            return p
        provenance['structural/'+aid] = sha(sq/'structural_qc_status.json')
        selection = stages.get(('selection', aid)); utr = stages.get(('utr', aid))
        if not selection or sha(next(p for k,a,p in request['stages'] if (k,a)==('selection',aid))) != ident['selection_sha256']:
            raise ValueError('Selection report differs from identified release')
        if selection.get('assembly_id') != aid or selection.get('reference_identity') != ref: raise ValueError('Selection report assembly/reference differs')
        if not utr or utr['reference_identity'] != ref: raise ValueError('Missing/mismatched UTR report')
        settings = by_assembly(preflight, aid)['settings']
        prediction = by_assembly(stages['braker','run'], aid)
        proteins = by_assembly(stages['protein','run'], aid)
        prepared = by_assembly(stages['rna','run'], aid)
        aligned = by_assembly(stages['rna_alignment','run'], aid)
        refstage = by_assembly(stages['reference_repeats','run'], aid)
        for state in (prediction,):
            if state.get('reference_identity') != ref: raise ValueError('Prediction/report reference differs')
        if prediction['status']!='ready' or refstage['status']!='ready': raise ValueError('Incomplete assembly report')
        has_rna=prediction['mode'] in ('rna','rna_protein');has_protein=prediction['mode'] in ('protein','rna_protein')
        if proteins['status']!=('ready' if has_protein else 'not_requested') or any(s['status']!=('ready' if has_rna else 'not_requested') for s in (prepared,aligned)):
            raise ValueError('Evidence mode/status mismatch')
        if refstage['reference']['reference_identity'] != ref: raise ValueError('Reference report differs')
        lengths = sorted((r['length'] for r in reference['contigs']), reverse=True); bases = sum(lengths)
        running = 0; n50 = None
        for n in lengths:
            running += n
            if running >= bases/2: n50 = n; break
        assembly_stats = refstage['reference'].get('stats', {})
        row = dict(assembly_id=aid, evidence_mode=prediction['mode'], genes=release['genes'], transcripts=release['coding_transcripts'],
            reference_bases=bases, retained_contigs=len(lengths), n50_bp=n50,
            repeat_union_bases=repeat.get('repeat_interval_union_bases'), functional_enabled=release['functional_enabled'])
        overview.append([aid,prediction['mode'],release['genes'],release['coding_transcripts'],f'[Release](release/{aid}/README.md)'])
        p = [f'## Assembly {ordinal}: {aid}', '', f'[Final deliverables](release/{aid}/README.md) · [Checksums](release/{aid}/checksums.sha256) · [Release manifest](release/{aid}/release_status.json)', '',
             '### Reference and evidence inputs', '', md_table(['Metric','Value'],[
                 ['Species',settings.get('species')],['NCBI taxid',settings.get('taxid')],['Evidence mode', prediction['mode']], ['Retained contigs',len(lengths)], ['Retained reference bases',bases],['Retained-reference N50 (bp)',n50],
                 ['Excluded contigs',assembly_stats.get('excluded_contigs')], ['Protein evidence records',proteins.get('protein_count')],
                 ['RNA readsets',len(prepared.get('rna_readsets',[]))]]), '',
             'N50 uses the retained reference contig lengths. Protein evidence records are database sequences, not predicted genes. Paired RNA read units are pairs; alignment fractions use primary individual reads.', '']
        if proteins.get('sources'):
            p += [md_table(['Protein source','Kind','Status','Protein records'],[[s['source_id'],s.get('request',{}).get('kind'),s['status'],s.get('protein_count')] for s in proteins['sources']]), '']
        ids = set(prepared.get('rna_readsets', []))
        read_units = {r['readset_id']:r for r in stages['rna','run'].get('readsets',[]) if r['readset_id'] in ids}
        mapping = [r for r in stages['rna_alignment','run'].get('readsets',[]) if r.get('assembly_id', r.get('strand',{}).get('assembly_id')) == aid]
        p += [md_table(['RNA readset','Selected units','Mapped primary reads','Resolved strand'],[[r['readset_id'],read_units.get(r['readset_id'],{}).get('read_units'),
            f"{r['mapped_fraction']:.1%}" if r.get('mapped_fraction') is not None else None,r.get('strand',{}).get('strandedness')] for r in mapping]) if mapping else 'RNA alignment: not requested.', '',
            '[Raw RNA MultiQC](rna/multiqc/raw/multiqc_report.html) · [Processed RNA MultiQC](rna/multiqc/filtered/multiqc_report.html) · [Alignment MultiQC](rna_alignment/multiqc/multiqc_report.html)' if ids else '', '',
            '### Prediction, UTR refinement and selection', '', md_table(['Model set','Genes','Coding transcripts'],[
                ['BRAKER predictions',prediction.get('genes'),prediction.get('coding_transcripts')],
                *[['Before selection' if r['model_set']=='before' else 'Final selected',r['genes'],r['transcripts']] for r in sqs['model_sets']]]), '',
            md_table(['UTR operation','Value'],[['Refinement skipped',utr['skipped']],['Changed transcripts',utr['changed_transcripts']],
                ['Rejected UTR proposals',len(utr.get('rejected_utr_candidates',{}))],['CDS sequences preserved',utr.get('cds_sequences_unchanged')],
                ['Protein sequences preserved',utr.get('protein_sequences_unchanged')]]), '',
            'Before/after below means accepted post-UTR models before selection versus the selected final-ID set.', '',
            f'![Model retention](qc/structural/{aid}/plots/retention.png)', '',
            'Selection rules and resolved thresholds:', '', md_table(['Setting','Value'],[[k,v] for k,v in sorted(selection['settings'].items())]), '',
            md_table(['Exclusion rule','Excluded transcripts triggering rule'],sorted(selection.get('removal_reason_counts',{}).items())) if selection.get('removal_reason_counts') else 'No transcript exclusion rules were triggered.', '',
            'One transcript can trigger multiple exclusion rules. Coding IDs and complete excluded-model decisions are available in the release ID maps.', '']
        p += [f'[Selection-decision figure](qc/structural/{aid}/plots/decisions.png)', '']
        completeness = tsv(qcfile('completeness_metrics.tsv'))
        for metric in ('busco_complete','busco_single','busco_duplicated','busco_fragmented','busco_missing','omark_consistent','omark_inconsistent','omark_contaminant','omark_unknown'):
            row[metric]=next((float(r['value']) for r in completeness if r['model_set']=='after' and r['metric']==metric),None)
        p += ['### Completeness and gene structure', '',
              md_table(['Model set','Metric','Value'],[['Before selection' if r['model_set']=='before' else 'Final selected',r['metric'].replace('busco_','BUSCO ').replace('omark_','OMArk ').replace('_',' ')+(' (orthogroups)' if r['metric']=='busco_total' else ' (%)'),int(float(r['value'])) if r['metric']=='busco_total' else float(r['value'])] for r in completeness]) if completeness else 'BUSCO/OMArk metrics: not requested.', '',
              'BUSCO and OMArk use one selected representative protein per coding gene. Percentages retain each tool’s native denominator; other isoforms are not independently classified.', '',
              f'![Completeness before and after selection](qc/structural/{aid}/plots/completeness.png)', '',
              f'[Gene/CDS structure distributions](qc/structural/{aid}/plots/structure.png) · [Intron/isoform distributions](qc/structural/{aid}/plots/complexity.png) · [AGAT and detailed structural QC](qc/structural/{aid}/README.md)', '']
        # Reuse accepted source-specific evidence; do not calculate a combined AED.
        assessment = tsv(payload('tables/transcript_assessment.tsv'))
        evidence_rows=[]
        for key,label in [('rna_aed','RNA AED'),('protein_aed','Protein AED'),('junction_coding_pass_fraction','Coding junction pass fraction'),('junction_utr_pass_fraction','UTR-associated junction pass fraction')]:
            values=[float(r[key]) for r in assessment if r.get(key) not in (None,'','.')] 
            if any(not math.isfinite(v) or not 0<=v<=1 for v in values): raise ValueError('Invalid evidence metric')
            evidence_rows.append([label,len(values),len(assessment)-len(values),statistics.median(values) if values else None])
            row[key+'_assessed']=len(values);row[key+'_median']=statistics.median(values) if values else None
        p += ['### RNA/protein evidence and splice-junction support', '', md_table(['Metric','Assessed transcripts','Unassessed / not applicable','Median'],evidence_rows), '',
              'RNA and protein AED are separate native InGenAnnot scores (lower means less discrepancy). Junction pass fractions are separate Portcullis-supported metrics. Missing evidence is not scored as AED=1 or zero junction support. Single-exon transcripts have no applicable intron fraction.', '',
              f'![Separate RNA and protein AED](qc/structural/{aid}/plots/aed.png)', '', f'![Coding and UTR splice support](qc/structural/{aid}/plots/junctions.png)', '',
              f'[Evidence agreement figure](qc/structural/{aid}/plots/evidence_agreement.png) · [Junction occurrences](qc/structural/{aid}/plots/junction_occurrences.png) · [Diagnostic flags](qc/structural/{aid}/plots/flags.png)', '',
              'Junction occurrence counts across transcripts can repeat shared genomic junctions. Native Portcullis reports retain distinct-junction and readset-level results.', '',
              '### Repeats and selected genome tracks', '', md_table(['Metric','Value'],[['Repeat interval-union bases',repeat.get('repeat_interval_union_bases')],
                  ['Repeat interval-union coverage',f"{repeat['repeat_interval_union_bases']/bases:.2%}" if repeat.get('repeat_interval_union_bases') is not None else None],
                  ['Lowercase A/C/G/T bases',repeat.get('lowercase_acgt_bases')],['Repeat alignment records',repeat.get('repeat_alignment_records')]]), '',
              'Coverage is an interval union. Repeat-class coverage can overlap and must not be summed. RepeatMasker RNA-labelled repeats are not added to tRNAscan/Rfam gene counts.', '',
              md_table(['Track','Status','Raw feature records','Selected feature records'],[[t,s['status'],s.get('raw_features'),s.get('features')] for t,s in track['tracks'].items()]), '']
        reasons=[]
        for t,s in track['tracks'].items():
            reasons += [[t,s.get('decision_unit'),k,n] for k,n in s.get('reasons',{}).items() if k not in ('passed','filter_disabled')]
        p += [md_table(['Track','Decision unit','Rejection reason','Count'],reasons) if reasons else 'No rejected feature/window decisions were recorded.', '',
              'Feature records can include exons and repeat fragments. Telomere decisions use candidate windows; selected regions merge adjacent qualifying windows. Centromeric satellite/HOR records are repeat predictions, not counts of functional centromeres.', '',
              f'[Raw calls, selected tracks, resolved filters and audits](genome_tracks/{aid}/README.md) · [Repeat preparation](reference_repeats/README.md)', '']
        present=[(t,s) for t,s in track['tracks'].items() if s['status']=='ready' and 'raw_features' in s]
        if present:
            fig,ax=plt.subplots(figsize=(8,3.8)); x=list(range(len(present)))
            ax.bar([i-.18 for i in x],[s['raw_features'] for t,s in present],.36,label='Raw feature records',color='#466A95')
            ax.bar([i+.18 for i in x],[s['features'] for t,s in present],.36,label='Selected feature records',color='#27856C')
            ax.set_xticks(x,[t for t,s in present]); ax.set_yscale('symlog',linthresh=1); ax.set_ylabel('Feature records (symlog scale)'); ax.set_title(aid+' — genome-track filtering');ax.legend()
            p += [save(fig,aid,'track_filtering'), '', 'Bars count feature records, not loci; telomere raw and selected regions have different merging boundaries.', '']
        if track['tracks']['telomere']['status']=='ready':
            telo=payload(f'tracks/{aid}.telomere.gff3'); candidates=[]
            for line in telo.read_text().splitlines():
                if not line.strip() or line.startswith('#'):continue
                c=line.split('\t'); a=dict(v.split('=',1) for v in c[8].split(';') if '=' in v)
                candidates.append([c[0],int(c[3]),int(c[4]),a.get('motif_density'),a.get('end_distance'),a.get('position')])
            p += [md_table(['Contig','Start','End','Motif density','End distance (bp)','Position'],candidates) if len(candidates)<=30 else f'{len(candidates):,} selected telomere regions; [full coordinates](release/{aid}/tracks/{aid}.telomere.gff3).', '']
            if candidates and all(c[3] is not None and c[4] is not None for c in candidates):
                fig,ax=plt.subplots(figsize=(7,3.8)); ax.scatter([int(c[4]) for c in candidates],[float(c[3])*100 for c in candidates],color='#27856C')
                ax.set(xlabel='Nearest contig-end distance (bp)',ylabel='Exact-motif coverage (%)',ylim=(0,100),title=aid+' — selected telomere regions')
                p += [save(fig,aid,'telomere_positions'), '', 'Each point is one selected merged region; overlapping points can represent several regions.', '']
        p += ['### Functional annotation', '']
        if release['functional_enabled']:
            functional=read(payload('provenance/functional_status.json'))
            p += [md_table(['Metric','Transcripts'],[['Named product',functional.get('named_transcripts')],['Hypothetical protein',functional.get('hypothetical_transcripts')]]), '',
                  md_table(['Tool','Hit','No hit','Ineligible','Disabled'],[[t,*[s.get(k,0) for k in ('hit','no_hit','ineligible','disabled')]] for t,s in functional.get('tools',{}).items()]), '',
                  'Coverage counts selected coding transcripts/proteins. Ineligible sequences remain in the release; disabled and no-hit results are distinct. OMArk classes apply only to assessed representatives.', '']
            tools=functional.get('tools',{})
            for tool,counts in tools.items():
                values=[counts.get(k,0) for k in ('hit','no_hit','ineligible','disabled')]
                if any(not isinstance(v,int) or v<0 for v in values) or sum(values)!=release['coding_transcripts']:
                    raise ValueError('Functional tool counts differ from selected transcripts: '+tool)
            row['named_transcripts']=functional.get('named_transcripts');row['hypothetical_transcripts']=functional.get('hypothetical_transcripts')
            if tools:
                fig,ax=plt.subplots(figsize=(8,4)); positions=list(range(len(tools))); left=[0]*len(tools)
                for category,color in [('hit','#27856C'),('no_hit','#466A95'),('ineligible','#B89046'),('disabled','#AAAAAA')]:
                    values=[s.get(category,0) for s in tools.values()]
                    ax.barh(positions,values,left=left,color=color,label=category.replace('_',' '));left=[a+b for a,b in zip(left,values)]
                ax.set_yticks(positions,list(tools));ax.set_xlabel('Selected coding transcripts');ax.set_title(aid+' — functional tool coverage');ax.xaxis.set_major_locator(MaxNLocator(integer=True));ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=4)
                p += [save(fig,aid,'functional_coverage'), '']
            # Carry forward the original report's vocabulary/product detail plots.
            # Count each term at most once per selected transcript (no enrichment).
            counters={k:Counter() for k in ('GO terms','InterPro domains','KEGG orthologs')}
            with payload('tables/functional_annotations.tsv').open() as handle:
                for record in csv.DictReader(handle,delimiter='\t'):
                    for label,keys,pattern in [('GO terms',('ipr_go','eggnog_GO'),r'GO:\d+'),
                                               ('InterPro domains',('ipr_accessions',),r'IPR\d+'),
                                               ('KEGG orthologs',('eggnog_KEGG_ko','kofam_KO'),r'K\d{5}')]:
                        counters[label].update(set(re.findall(pattern,' '.join(record.get(k,'') for k in keys))))
            products=Counter()
            with payload('tables/gene_product_summary.tsv').open() as handle:
                for record in csv.DictReader(handle,delimiter='\t'):products[unquote(record['product'])]+=1
            counters['Gene products']=products
            links=[]
            for label,counter in counters.items():
                if not counter:continue
                top=sorted(counter.items(),key=lambda kv:(-kv[1],kv[0]))[:20];stem='top_'+label.lower().replace(' ','_')
                labels=[textwrap.fill(k,50) for k,v in reversed(top)]
                fig,ax=plt.subplots(figsize=(9, max(3,sum(k.count('\n')+1 for k in labels)*.25+1.1)))
                ax.barh(labels,[v for k,v in reversed(top)],color='#466A95');ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax.set_xlabel('Selected coding genes' if label=='Gene products' else 'Selected coding transcripts with term')
                ax.set_title(aid+' — '+label)
                save(fig,aid,stem)
                links.append(f'[{label}](reports/annotation_report_plots/{aid}.{stem}.png)')
            if links:p += ['Detailed functional figures: '+' · '.join(links), '',
                          'Vocabulary plots count each term once per selected transcript; product plots count selected genes. They show annotation frequency, not enrichment or biological interpretation.', '']
            p += [f'[Functional searches, products and detailed results](functional_annotation/{aid}/README.md)', '']
        else:
            row['named_transcripts']=None;row['hypothetical_transcripts']=None
            p += ['Functional annotation: disabled.', '']
        p += ['### Methods, settings and source reports', '', md_table(['Stage','Status'],[
            ['Protein preparation',proteins['status']],['RNA preparation',prepared['status']],['RNA alignment',aligned['status']],
            ['BRAKER prediction',prediction['status']],['UTR refinement','skipped' if utr['skipped'] else 'ready'],
            ['Model selection',selection['status']],['Final structural QC',sqs['status']],
            ['Native AED','enabled' if request.get('configuration',{}).get('run_aed') is True else 'disabled' if request.get('configuration',{}).get('run_aed') is False else 'not recorded'],
            ['Junction QC','enabled' if request.get('configuration',{}).get('run_junction_qc') is True else 'disabled' if request.get('configuration',{}).get('run_junction_qc') is False else 'not recorded'],
            ['Release packaging',release['status']]]), '',
            'Configured/resolved settings below describe the requested execution. They do not establish the release of an already-populated shared database; `CURRENT` is an acquisition setting. Native runtime/version records and source manifests remain linked with each stage.', '',
            md_table(['Per-assembly setting','Value'],[[k,json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v] for k,v in sorted(settings.items())]), '',
            md_table(['Genome-track filter','Resolved policy'],[[t,json.dumps(s.get('filter_policy',{}),sort_keys=True)] for t,s in track['tracks'].items() if s['status']=='ready']) if any(s['status']=='ready' for s in track['tracks'].values()) else 'All optional genome tracks were disabled.', '',
            f'[Prediction and native runtime logs](structural_annotation/README.md) · [UTR changes and checks](structural_annotation/utr/{aid}/README.md) · [Full structural QC and figures](qc/structural/{aid}/README.md)', '',
            '[Preflight and taxonomy provenance](preflight/README.md) · [Protein evidence sources](protein_evidence/README.md) · [RNA preparation](rna/README.md) · [RNA alignment](rna_alignment/README.md)', '',
            f'Reference identity: `{ref}`.', '']
        # Inline structural figures must exist and match the accepted QC manifest.
        for name in ('retention','completeness','aed','junctions','structure','complexity','evidence_agreement','junction_occurrences','flags','decisions'):
            qcfile(f'plots/{name}.png')
        parts += p; metrics.append(row)
    header=['# Genome Annotation Report', '', '[Final release index](release/README.md)', '',
        md_table(['Assembly','Evidence mode','Coding genes','Coding transcripts','Deliverables'],overview), '',
        '## Contents', '', *[f'- [Assembly {i}: {escape(a)}](#assembly-{i}-{re.sub(r"[^a-z0-9_-]", "", a.lower())})' for i,a in enumerate(expected,1)], '',
        'Results are reported per assembly. Missing/not-requested evidence remains distinct from zero or no hit. This report does not calculate a composite confidence score, apply additional filtering or provide tuning recommendations.', '']
    text='\n'.join(header+parts+['## Run configuration', '', md_table(['Configured parameter','Value'],[[k,json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v] for k,v in sorted(request.get('configuration',{}).items())]), '',
        '[Machine-readable report status](reports/annotation_report_status.json) · [Assembly metrics](reports/annotation_metrics.tsv) · [Source stage manifests](reports/annotation_report_data/)', ''])
    (out/'annotation_report.md').write_text(text,encoding='utf-8')
    with (out/'annotation_metrics.tsv').open('w',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=list(metrics[0]) if metrics else ['assembly_id'],delimiter='\t');writer.writeheader();writer.writerows(metrics)
    output_paths=[out/'annotation_report.md',out/'annotation_metrics.tsv']+[p for folder in (assets,data) for p in folder.rglob('*') if p.is_file()]
    outputs={p.relative_to(out).as_posix():sha(p) for p in output_paths}
    status=dict(schema_version=1,stage='annotation_report',status='ready',assemblies=expected,source_sha256=provenance,outputs=outputs,
        figure_stems=figures,report_policy='One Markdown report; assembly subsections; key figures inline; results and methods only; no pipeline HTML or tuning advice')
    dump(out/'annotation_report_status.json',status)
    return status


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--request',required=True);parser.add_argument('--outdir',required=True)
    args=parser.parse_args();build(json.loads(base64.b64decode(args.request)),args.outdir)


if __name__=='__main__':main()
