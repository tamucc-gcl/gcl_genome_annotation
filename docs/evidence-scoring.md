# Native RNA and protein AED

Production uses InGenAnnot 0.0.20 directly to assess accepted post-UTR models.
It reports **`rna_aed`** and **`protein_aed`**, selected evidence IDs and the
RNA penalty indicator. Lower values indicate closer agreement with the evidence
chosen by that method; they are not probabilities of gene correctness.
Missing evidence is unavailable (`.`), never an invented AED of 1 or 0.
There is no production combined AED, duplicate alias or custom legacy AED.

AED is a disagreement measure using sensitivity and specificity:
`AED = 1 - (sensitivity + specificity) / 2`. Its meaning depends on the
evidence preparation, coordinate scope, chain selection and native penalties,
not just that arithmetic. InGenAnnot is the selected published implementation;
the pipeline does not claim to reproduce MAKER's evidence clustering or scores.

## Evidence preparation

Protein evidence uses the existing miniprot alignment and complete supplied
prediction protein sources. RNA uses retained BRAKER ETP/StringTie merged GTF
when present and compatible; otherwise the existing StringTie assembly/merge
modules use the separate accepted BAMs. `aed_reassemble_rna=true` forces that
fallback. RNA evidence is shared with UTR refinement. Protein alignment depends
on reference/proteins, not UTR settings, and can overlap prediction.

`ingenannot_rna_cds_only=true` and `ingenannot_rna_penalty=0.25` are the default
profile. Protein AED uses native CDS scope. RNA AED therefore assesses coding
support rather than complete UTR boundary quality. Native evidence selection and
overflow penalty remain explicit. Large protein sets can produce optimistic
best-match support; neither a broad nor a narrow score distribution alone
establishes biological quality. The default penalty is not a calibrated model
selection threshold.

## Integrity and outputs

The native output is used for scores, not as permission to change model
structures. The wrapper validates model/reference/evidence identities and coding
structure, distinguishing known native serialization/bounds normalization from
real changes. Production outputs are under
`structural_annotation/aed/<assembly>/`: scores, status, evidence provenance and
Markdown. Structural assessment joins them to the exact transcript population;
final selection re-keys retained rows without rerunning AED.

`run_aed=false` disables assessment. AED filtering is separately opt-in using
`filter_aed_source=rna|protein`, `filter_max_aed`, and the chosen action. Current
selection is per transcript with gene rescue by a surviving isoform; absent
evidence for a mode is explicitly skipped. See [selection](model-selection.md).
Portcullis junction QC is complementary and does not alter the native AED formula.

The custom scorer and unchanged-model comparison entries have been retired.
Production uses `modules/ingenannotAed.nf`, `ingenannot_inputs.py` and
`ingenannot_aed.py`; common GFF/GTF/interval validation lives in
`annotation_models.py`. No custom AED arithmetic or comparison output is retained
in production. Old benchmark/custom-score directories can persist from earlier
runs; they are not current results. See [consolidation record](cleanup.md).

[Native AED documentation](https://bioger.pages-forge.inrae.fr/ingenannot/tools/AED.html),
[InGenAnnot study](https://pubmed.ncbi.nlm.nih.gov/40971313/),
[UTR integration](ingenannot.md), [tests](testing.md).
