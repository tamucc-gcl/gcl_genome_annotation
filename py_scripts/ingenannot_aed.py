"""Validate unchanged native InGenAnnot assessment and read its RNA/protein scores."""
import json
from pathlib import Path
from annotation_models import features, models
from prepare_protein_evidence import checksum

def assess(original, annotated, identity, inputs):
    identity = json.loads(Path(identity).read_text())
    accepted = json.loads(Path(inputs).read_text())
    if accepted.get('status') != 'ready' or identity['reference_identity'] != accepted['reference_identity'] or checksum(original) != accepted['model_sha256']:
        raise ValueError('Native AED reference/models differ from validated inputs')
    lengths = {c['id']:c['length'] for c in identity['contigs']}
    before, after = models(original,lengths), models(annotated,lengths)
    if set(before) != set(after):
        raise ValueError('InGenAnnot changed transcript identities during assessment')
    normalized_bounds = []
    repeated_exons = {}
    for tid,a in before.items():
        b = after[tid]
        # Native export can repeat identical exon records for shared isoforms.
        # Compare exact exon boundaries, not their merged coverage: splitting or
        # joining an exon must still fail. CDS records and phases remain strict.
        if (any(a[k]!=b[k] for k in ('gene_id','seq','strand'))
                or set(a['exons']) != set(b['exons'])
                or sorted(a['cds']) != sorted(b['cds'])):
            raise ValueError(f'InGenAnnot changed model structure during assessment: {tid}')
        if sorted(a['exons']) != sorted(b['exons']):
            repeated_exons[tid] = dict(original_records=len(a['exons']),native_records=len(b['exons']),
                                       unique_exons=len(set(a['exons'])))
        original_bounds = (a['start'],a['end'])
        native_bounds = (b['start'],b['end'])
        if original_bounds != native_bounds:
            # InGenAnnot 0.0.20 Utils fixes transcript bounds to the exon envelope.
            # BRAKER may include a separately recorded stop codon in its mRNA span.
            exon_bounds = (min(s for s,e in a['exons']),max(e for s,e in a['exons']))
            if native_bounds != exon_bounds:
                raise ValueError(f'InGenAnnot changed model structure during assessment: {tid}: unexpected transcript bounds')
            normalized_bounds.append(tid)
    def coding_phases(path):
        result = {tid:[] for tid in before}
        for row in features(path,lengths):
            if row['kind'] == 'CDS':
                for tid in row['attrs']['Parent'].split(','):
                    result[tid].append((row['start'],row['end'],row['phase']))
        return {tid:sorted(parts) for tid,parts in result.items()}
    if coding_phases(original) != coding_phases(annotated):
        raise ValueError('InGenAnnot changed model structure during assessment: coding phases differ')
    scores = {}
    for row in features(annotated,lengths):
        if row['kind'] not in ('mRNA','transcript'): continue
        attrs = row['attrs']
        tid = attrs['ID']
        entry = {}
        for kind,field,evidence in [('rna','aed_ev_tr','ev_tr'),('protein','aed_ev_pr','ev_pr')]:
            value = float(attrs[field])
            if not 0 <= value <= 1:
                raise ValueError('Invalid native AED value')
            entry[kind] = value if accepted['evidence'][kind]['available'] else None
            entry[kind+'_evidence_id'] = attrs.get(evidence)
        entry['rna_penalty'] = attrs.get('ev_tr_penalty')
        scores[tid] = entry
    accepted['native_exon_record_normalizations'] = repeated_exons
    return before,scores,normalized_bounds,accepted
