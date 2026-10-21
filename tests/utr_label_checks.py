#!/usr/bin/env python3
"""UTR boundary-label regression helpers using native exon/CDS geometry."""
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from annotation_models import features, merged, overlap_length
from utr_refinement import subtract


def annotations(path,lengths):
    result={}
    def entry(tid):
        return result.setdefault(tid,dict(exons=[],cds=[],stops=[],five_prime_UTR=[],three_prime_UTR=[],feature_contexts=set()))
    for row in features(path,lengths):
        kind=row['kind'];attrs=row['attrs']
        if kind in ('mRNA','transcript'):
            entry(attrs['ID']).update(seq=row['seq'],strand=row['strand'],attributes=attrs)
        elif kind in ('exon','CDS','stop_codon','five_prime_UTR','three_prime_UTR'):
            key={'exon':'exons','CDS':'cds','stop_codon':'stops'}.get(kind,kind)
            for tid in attrs['Parent'].split(','):
                e=entry(tid);e[key].append((row['start'],row['end']))
                e['feature_contexts'].add((row['seq'],row['strand']))
    return result


def inspect(original,native,identity,output):
    reference=json.loads(Path(identity).read_text()); lengths={c['id']:c['length'] for c in reference['contigs']}
    print('Reading original coding structures...',flush=True)
    old=annotations(original,lengths)
    print('Reading completed native UTR output...',flush=True)
    new=annotations(native,lengths); issues=[]; eligible_issues=0
    for tid,b in new.items():
        a=old[tid]
        if not b['cds']: continue
        lower=min(s for s,e in b['cds']);upper=max(e for s,e in b['cds'])
        low=[(s,min(e,lower)) for s,e in b['exons'] if s<lower]
        high=[(max(s,upper),e) for s,e in b['exons'] if e>upper]
        expected={'five_prime_UTR':low if b['strand']=='+' else high,'three_prime_UTR':high if b['strand']=='+' else low}
        rejection=[]
        if b['feature_contexts']!={(b['seq'],b['strand'])} or (a['seq'],a['strand'])!=(b['seq'],b['strand']):rejection.append('contig_or_strand')
        if sorted(a['exons'])!=sorted(b['exons']):
            if not b['attributes'].get('utr_refine_evidence'):rejection.append('missing_provenance')
            if overlap_length(a['exons'],b['exons'])!=sum(e-s for s,e in merged(a['exons'])):rejection.append('shortening')
            if subtract(a['exons'],a['cds']):rejection.append('existing_utr_change')
            if overlap_length(subtract(b['exons'],a['exons']),[(lower,upper)]):rejection.append('internal_coding_change')
        for kind,blocks in expected.items():
            want=subtract(blocks,a['stops']);actual=subtract(b[kind],a['stops'])
            if want==actual:continue
            if not rejection:eligible_issues+=1
            issues.append(dict(transcript_id=tid,kind=kind,seq=b['seq'],strand=b['strand'],
                evidence_id=b['attributes'].get('utr_refine_evidence'),candidate_rejection_reasons=rejection,
                expected=want,native=actual,missing_label_blocks=subtract(want,actual),extra_label_blocks=subtract(actual,want),
                original_exons=a['exons'],native_exons=b['exons'],cds=b['cds'],stop_codons=a['stops'],
                boundary_equalities=dict(cds_start_is_exon_last_base=any(e==lower+1 for s,e in b['exons']),
                    cds_end_is_exon_first_base=any(s==upper-1 for s,e in b['exons']))))
    result=dict(stage='utr_label_diagnostic',coordinates='0-based half-open',
        scope='Read-only comparison of completed native labels with exon/CDS geometry; no refinement, ranking or sequence extraction rerun',
        transcripts=len(new),mismatched_labels=len(issues),mismatched_labels_on_eligible_candidates=eligible_issues,issues=issues)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(f"Found {len(issues)} label mismatches, {eligible_issues} on otherwise eligible candidates. Report: {output}",flush=True)
    return result


