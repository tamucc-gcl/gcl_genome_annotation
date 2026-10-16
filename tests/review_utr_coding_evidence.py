#!/usr/bin/env python3
"""Compare coding/protein support at UTR/CDS overlaps from completed outputs."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from evidence_aed import models, merged, nucleotide_aed, overlap_length
from prepare_protein_evidence import checksum, records
from review_ingenannot_evidence import table, selected_proteins


def compare_chain(chain,model):
    if (chain['seq'],chain['strand'])!=(model['seq'],model['strand']):
        return dict(overlapping_cds_bp=0,model_cds_fraction=0)
    common=overlap_length(chain['blocks'],model['cds'])
    length=sum(e-s for s,e in merged(model['cds']))
    return dict(overlapping_cds_bp=common,model_cds_fraction=common/length if length else None)


def review(root,aid,review_path,outdir):
    root,outdir=Path(root),Path(outdir);base=root/'structural_annotation'
    candidates=json.loads(Path(review_path).read_text())
    benchmark=base/'ingenannot_benchmark'/aid
    status=json.loads((benchmark/'benchmark_status.json').read_text())
    original=base/'braker'/aid/f'{aid}.gff3'
    identity=json.loads((root/'filtered_assemblies'/aid/f'{aid}.contigs.json').read_text())
    if status['status']!='ready' or status['model_sha256']!=checksum(original) or status['reference_identity']!=identity['reference_identity']:
        raise ValueError('Coding review requires scores from the accepted original models/reference')
    parsed=models(original,{c['id']:c['length'] for c in identity['contigs']})
    pairs=[];genes=set()
    for case in candidates['selected_cases']:
        for neighbor in case['other_gene_coding_spans_overlapped']:
            if not neighbor['extension_overlap_cds_bp']:continue
            tid=case['transcript_id'];model=parsed[tid]
            if (model['seq'],model['strand'],sorted(model['cds']))!=(case['original']['seq'],case['original']['strand'],sorted(map(tuple,case['original']['cds']))):
                raise ValueError('Candidate review differs from original coding models')
            gid=model['gene_id'];other=neighbor['gene_id']
            genes.update((gid,other))
            pairs.append(dict(transcript_id=tid,gene_id=gid,other_gene_id=other,
                same_strand=model['strand']==neighbor['strand'],utr_overlap_cds_bp=neighbor['extension_overlap_cds_bp']))
    chosen={tid:m for tid,m in parsed.items() if m['gene_id'] in genes}
    scores=table(benchmark/'benchmark_scores.tsv')
    wanted={scores[t]['protein_evidence_id'] for t in chosen}-{'','.','None'}
    print(f'Reading selected protein chains for {len(chosen)} coding models. One streaming scan; no alignment...',flush=True)
    chains=selected_proteins(base/'evidence_aed'/aid/'protein_alignment'/f'{aid}.miniprot.gff',wanted,include_attributes=True)
    # Retain lengths, not full sequences, for a compact report.
    proteins={}
    with (base/'braker'/aid/f'{aid}_proteins.fa').open() as handle:
        for header,sequence in records(handle):
            tid=header.split()[0]
            if tid in chosen:proteins[tid]=dict(length_aa=len(sequence),sequence_ends_with_stop=sequence.endswith('*'))
    if set(proteins)!=set(chosen):raise ValueError('Selected predicted protein IDs missing')
    details=[]
    for tid,model in sorted(chosen.items()):
        score=scores[tid];chain=chains.get(score['protein_evidence_id'])
        measured=None
        if chain:
            if (chain['seq'],chain['strand'])!=(model['seq'],model['strand']):raise ValueError('Selected protein chain has incompatible contig/strand')
            measured=nucleotide_aed(model['cds'],chain['blocks'])['aed']
            if abs(measured-float(score['ingenannot_protein_aed']))>.000051:
                raise ValueError('Selected protein chain differs from accepted protein AED')
        details.append(dict(transcript_id=tid,gene_id=model['gene_id'],contig=model['seq'],strand=model['strand'],
            cds=merged(model['cds']),cds_bp=sum(e-s for s,e in merged(model['cds'])),predicted_protein=proteins[tid],
            native_rna_aed=score['ingenannot_rna_aed'],native_protein_aed=score['ingenannot_protein_aed'],
            recomputed_protein_aed=measured,selected_protein_chain=chain,
            selected_chain_overlap_other_coding_models={other:compare_chain(chain,m) for other,m in chosen.items()
                if chain and m['gene_id']!=model['gene_id'] and m['seq']==model['seq']}))
    result=dict(stage='utr_coding_evidence_review',coordinates='0-based half-open',
        model_sha256=checksum(original),reference_identity=identity['reference_identity'],pairs=pairs,models=details,
        scope='Selected coding structures and existing protein alignments; not independent validation or an automatic UTR rejection rule')
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'coding_review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'Coding evidence report: {outdir}/coding_review.json',flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs-root',required=True);p.add_argument('--assembly-id',required=True)
    p.add_argument('--review',required=True);p.add_argument('--outdir',required=True)
    a=p.parse_args();review(a.inputs_root,a.assembly_id,a.review,a.outdir)
