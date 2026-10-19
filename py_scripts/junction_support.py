#!/usr/bin/env python3
"""Match native Portcullis junctions to accepted exon boundaries. Report only."""
import argparse
import csv
import json
from pathlib import Path

from evidence_aed import models
from portcullis_inputs import decode, dump, table
from prepare_protein_evidence import checksum


def write_table(path,columns,rows):
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({k:'.' if row[k] is None else row[k] for k in columns})


def report(gff,validation,identity,request,folders,outdir):
    request=decode(request);state=json.loads(Path(validation).read_text());reference=json.loads(Path(identity).read_text())
    if (state['status']!='ready' or state.get('stage')!='final_annotation'
            or state['gff3_sha256']!=checksum(gff) or state['reference_identity']!=reference['reference_identity']):
        raise ValueError('Junction QC requires unchanged accepted final annotation and reference')
    aid=request['assembly_id'];expected={r['readset_id']:r for r in request['readsets']}
    if len(expected)!=len(request['readsets']):raise ValueError('Duplicate requested readsets')
    lengths={c['id']:c['length'] for c in reference['contigs']}
    annotation=models(gff,lengths)
    libraries={};provenance=[]
    for folder in folders:
        folder=Path(folder);native=json.loads((folder/'junctions.json').read_text())
        rid=native['readset']['readset_id']
        if (rid not in expected or rid in libraries or native['assembly_id']!=aid
                or native['status']!='ready' or native['reference_identity']!=reference['reference_identity']
                or native['readset']!=expected[rid]):
            raise ValueError('Junction readset/reference metadata mismatch')
        for name,digest in native['table_sha256'].items():
            if checksum(folder/name)!=digest:raise ValueError('Native junction table changed after validation')
        raw=table(folder/'all.junctions.tab',lengths)
        passed=table(folder/'filtered.pass.junctions.tab',lengths)
        libraries[rid]=(raw,passed)
        provenance.append(native)
    if set(libraries)!=set(expected):raise ValueError('Missing/unexpected native readsets for junction QC')
    assessed=bool(expected);junction_rows=[];transcripts=[]
    for tid,model in sorted(annotation.items()):
        exons=sorted(set(model['exons']));coding=sorted(set(model['cds']))
        if any(left[1]>=right[0] for left,right in zip(exons,exons[1:])):
            raise ValueError(f'Overlapping/adjacent exons do not define a valid intron: {tid}')
        details=[]
        for left,right in zip(exons,exons[1:]):
            start,end=left[1],right[0]
            # Coding junctions connect coding sequence on BOTH exon sides.
            # Everything else is a UTR-associated junction (possibly mixed).
            is_coding=(any(s<left[1] and e==left[1] for s,e in coding)
                       and any(s==right[0] and e>right[0] for s,e in coding))
            key=(model['seq'],start,end,model['strand']);ambiguous=key[:3]+('.',)
            observed=[];passing=[];failed=[];unknown=[];alignments=0
            for rid,(raw,good) in sorted(libraries.items()):
                if key in raw:
                    observed.append(rid);alignments+=raw[key]['raw_alignments']
                    (passing if key in good else failed).append(rid)
                if ambiguous in raw:unknown.append(rid)
            label=('not_assessed' if not assessed else 'pass' if passing else 'fail' if observed
                   else 'unresolved_strand' if unknown else 'not_observed')
            row=dict(gene_id=model['gene_id'],transcript_id=tid,contig=model['seq'],strand=model['strand'],
                intron_start=start+1,intron_end=end,region='coding' if is_coding else 'utr_associated',
                status=label,passing_readsets=','.join(passing) or '.',failed_readsets=','.join(failed) or '.',
                unresolved_strand_readsets=','.join(unknown) or '.',raw_alignments=alignments if assessed else None)
            junction_rows.append(row);details.append(row)
        tx=dict(gene_id=model['gene_id'],transcript_id=tid,
            status='no_rna' if not assessed else 'single_exon' if not details else 'assessed')
        for prefix,selected in [('all',details),('coding',[r for r in details if r['region']=='coding']),
                                ('utr',[r for r in details if r['region']=='utr_associated'])]:
            tx[prefix+'_junctions']=len(selected)
            for label in ('pass','fail','not_observed','unresolved_strand'):
                tx[prefix+'_'+label]=sum(r['status']==label for r in selected) if assessed else None
            tx[prefix+'_pass_fraction']=sum(r['status']=='pass' for r in selected)/len(selected) if assessed and selected else None
        transcripts.append(tx)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    columns=['gene_id','transcript_id','contig','strand','intron_start','intron_end','region','status',
             'passing_readsets','failed_readsets','unresolved_strand_readsets','raw_alignments']
    write_table(out/'junction_support.tsv',columns,junction_rows)
    write_table(out/'transcript_junction_support.tsv',list(transcripts[0]),transcripts)
    summary=dict(status='ready',stage='junction_qc',assembly_id=aid,assessed=assessed,
        reason=None if assessed else 'No supplied RNA evidence',model_sha256=checksum(gff),
        reference_identity=reference['reference_identity'],transcripts=len(transcripts),readsets=sorted(provenance,key=lambda r:r['readset']['readset_id']),
        transcript_junctions=len(junction_rows),unique_model_junctions=len({(r['contig'],r['strand'],r['intron_start'],r['intron_end']) for r in junction_rows}),
        matching='Exact intron boundaries and resolved strand; pass in any readset; unresolved strands reported separately',
        coordinates='1-based inclusive intron coordinates',
        scope='Read-only junction QC; no BAM/model filtering, annotation rescue or AED modification')
    summary['counts']={label:sum(r['status']==label for r in junction_rows)
                       for label in ('pass','fail','not_observed','unresolved_strand','not_assessed')}
    dump(out/'junction_status.json',summary)
    (out/'README.md').write_text('# RNA splice-junction support\n\nStatus: **ready**\n\n'
        f"RNA assessed: {assessed}. Transcripts: {len(transcripts)}. Model junction occurrences: {len(junction_rows)}.\n\n"
        '[Transcript summary](transcript_junction_support.tsv) | [Junction details](junction_support.tsv) | [Provenance](junction_status.json)\n\n'
        'Portcullis runs separately for each readset. A junction passes if at least one readset passes it. '
        'Fail means observed on the matching strand but passing in no readset; not_observed means no matching junction. '
        'Unknown-strand observations are reported separately and never counted as passing strand-specific support. '
        'Raw counts are spliced alignments, not independent fragments.\n\n'
        'Coding junctions connect CDS on both sides; other junctions are UTR-associated. '
        'Single-exon transcripts and categories with no junctions have no support fraction. '
        'RNA-free assemblies are not assessed. Missing support does not prove an incorrect model. '
        'Small datasets may trigger Portcullis native rule-based fallback; inspect readset logs and provenance.\n\n'
        'No predicted annotation is supplied as a rescue reference. BAMs, gene models and AED remain unchanged. '
        'Individual supported junctions do not establish that the entire transcript chain occurs in one RNA molecule.\n\n'
        '[Final annotation](../../../structural_annotation/final/'+aid+'/README.md)\n')
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('gff','validation','identity','request','outdir'):parser.add_argument('--'+name,required=True)
    parser.add_argument('--readset',action='append',default=[])
    args=parser.parse_args()
    report(args.gff,args.validation,args.identity,args.request,args.readset,args.outdir)
