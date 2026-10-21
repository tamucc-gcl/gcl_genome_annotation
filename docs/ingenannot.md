# InGenAnnot UTR refinement and integration

Production uses pinned InGenAnnot 0.0.20 for optional UTR refinement and native
RNA/protein AED. The earlier custom-score benchmark is complete as a method
selection exercise; it is not run on every annotation. See
[AED semantics](evidence-scoring.md) for the implemented scoring profile.

The connected order is:

`BRAKER -> optional UTR refinement -> validated unfiltered models -> unfiltered QC/AED/junction assessment -> selection -> final IDs and selected QC -> function -> release/report`

## RNA and native ranking

Retained ETP StringTie merged transcripts are reused when available; otherwise
StringTie 2.2.1 assembles/merges separate accepted readset BAMs. The same RNA
evidence serves UTR refinement and AED. There is no BAM merge. Native
`clusterize -g`, `isoform_ranking` with paired/strand FOF metadata and
`utr_refine --utr_mode rank --onlynew` choose compatible structures. The pipeline
does not choose the longest RNA transcript unconditionally.

`utr_refine=true` is the default; RNA-free and explicitly disabled assemblies
pass their original models through with a skipped status. The pinned ranking
interface supports reverse/unstranded libraries; forward/auto metadata is
explicitly blocked until a validated adapter exists. Use `--utr_refine false`
for known forward-stranded inputs; mapping and native RNA AED can still run.
Default ranking requests are 16 CPUs/64 GB/96 h, with `utr_sj_threshold=0.05`,
`utr_cov_threshold=0.05`, `utr_alt_threshold=0.1`. These are configurable native
settings. `environments/ingenannot.yaml` is active production configuration.

## Safeguards and native adaptation

- Preserve original BRAKER IDs/parents, CDS coordinates/phases and separately
  annotated start/stop codons. gffread 0.12.7 checks baseline/refined CDS and
  protein sequence invariance; accepted UTR changes must not change coding content.
- Optional proposals with wrong contig/strand, shortening, changes to existing
  UTRs, missing native evidence provenance or added exon bases within the coding
  span restore the original transcript and record a rejection. A native coding
  change is checked before restoration and remains a validation failure.
- Duplicate exon records with exactly the original coordinate set are normalized
  to the original representation and audited, not miscounted as biological
  extensions. Other segmentation/base changes are not accepted by that rule.
- Remove separately annotated stop-codon bases from native UTR labels. Complete
  missing labels only for the known pinned-tool coding-boundary equality cases,
  using native-selected exons and unchanged CDS. Unexpected labels still fail.
- `utr_warnings.tsv` records accepted UTR additions overlapping other genes' CDS,
  including strand/overlap context. These warnings do not filter models or prove
  coding assignment. There is no species-specific UTR length exclusion.

## Outputs

`structural_annotation/utr/<assembly>/` contains validated refined GFF3,
`utr_changes.tsv`, `utr_warnings.tsv`, `utr_status.json`, native candidates/ranking
under `native/`, rejection/normalization audits and unfiltered extracted
`transcripts.fa`. This stage's transcript FASTA is not the final selected-ID
release FASTA. Main publishes accepted unfiltered models at
`structural_annotation/final/<assembly>/` and later exports selected final-ID
transcripts through [release packaging](annotation-release.md).

The original BRAKER protein/CDS files are retained as canonical coding sequences.
Validation/report helpers are separated from expensive ranking and native AED,
so unchanged scientific inputs can reuse those tasks after reporting changes.
Unfiltered and selected QC use the common AGAT/BUSCO/OMArk workflow.

## Tests and focused entries

```bash
python gcl_genome_annotation/tests/validate_utr_cluster.py \
  --outdir testing/utr_contracts_01
python gcl_genome_annotation/tests/validate_postprocessing_cluster.py \
  --outdir testing/postprocessing_contracts_01
```

Native cluster/rank/refine, coding invariance, both strands of model coordinates,
disabled/no-RNA behavior and caching passed synthetic acceptance. Model strand
is distinct from the library-strand restriction above. The connected four-mode
postprocessing fixture and real UTR/AED pilot were reviewed separately. Targeted
RNA support checks do not prove every transcript boundary/coding assignment.

`utr_refinement.nf` reads accepted published inputs via `utr_inputs_root` and
stops after validation and standalone UTR AGAT. Main uses native channels and
does not require that invocation. The completed benchmark entry was removed;
its production preparation/run/validation code is consolidated in
`modules/ingenannotAed.nf` and `py_scripts/ingenannot_aed.py`.
See [cleanup](cleanup.md) for retired names and cache implications.

[Native UTR workflow](https://bioger.pages-forge.inrae.fr/ingenannot/usecases/add_utrs.html),
[UTR command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/UTRRefine.html),
[study](https://pubmed.ncbi.nlm.nih.gov/40971313/).
