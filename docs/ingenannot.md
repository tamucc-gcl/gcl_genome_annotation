# InGenAnnot benchmark and UTR integration

The user selected a peer-reviewed existing evaluator rather than requiring exact
MAKER equivalence. InGenAnnot 0.0.20 is the candidate. Its native AED command
provides RNA and protein scores; combined AED is omitted. The current custom
method is retained as a benchmark baseline, not silently replaced.

## First stage: unchanged-model benchmark

### Selected-case coordinate review

After a completed real benchmark, run this from the project directory:

```bash
python gcl_genome_annotation/tests/review_ingenannot_evidence.py \
  --inputs-root testing/genome_annotation \
  --assembly-id Sfl-n_primary \
  --outdir testing/ingenannot_review
```

This selects up to two distinct genes per category: changed RNA choice, RNA
penalty, weak RNA/protein support with selected evidence, native score 1 despite
previous overlap, and strong agreement controls. It checks annotation/RNA
identities, recomputes the CDS-only RNA comparison and overflow conditions for
all overlapping same-strand assembled RNA candidates at each selected locus,
and compares the best penalized result with the native score. Selected protein
chains are extracted in one streaming scan of the existing miniprot GFF and
their structural scores checked. No alignments or native scoring are rerun.

`review.json` contains coordinate blocks, selected/previous RNA, top candidates,
CDS bases outside RNA, RNA extra exonic bases inside the CDS span, matching
assembled exon-boundary introns, and score agreement checks. Coordinates are
zero-based half-open. These intron comparisons do not measure independent read
junction support. A disagreement is recorded for inspection rather than hidden.
The utility requires the retained ETP RNA assembly used by this real benchmark.

Summary metrics use `score_equals_one`, `no_selected_evidence`, and
`score_one_with_previous_overlap`; native score 1 is not labeled zero overlap.

Input validation lives in `ingenannot_inputs.py`; preparation stages only that
helper and its checksum dependency. Reporting stages `ingenannot_benchmark.py`
and its own dependencies. Report edits therefore do not invalidate preparation
or native scoring. Installing this split changes the preparation task once;
existing completed native scores can still be reported directly without rerunning
assessment. The cluster contract edits a private copy of the report helper and
requires preparation/scoring cache hits and a new report task.

InGenAnnot 0.0.20 recalculates transcript bounds from exon bounds on import.
BRAKER can include a separately recorded stop codon outside its exon/CDS span:
the real `g1.t1` mRNA starts at 4259, but its first exon/CDS starts at 4262.
The report accepts only unchanged bounds or the exact unchanged exon envelope,
and records affected transcript IDs in `transcript_bounds_normalized_to_exons`.
Transcript IDs, gene parents, reference, strand, exact exon/CDS blocks and CDS
phases must remain unchanged. The original BRAKER annotation is retained;
the native GFF is an assessment export, not a replacement annotation.
This reporting correction does not change native scoring and can reuse its
completed task with `-resume`.

The first real run demonstrated slow full-database scaling: after about 21 hours
with four workers, RNA had finished 216/216 batches and protein 61/216. Batch
costs vary; these counts are not a reliable precise ETA. The original 24-hour
task limit was insufficient. RUN_INGENANNOT_AED now defaults to 64 CPUs,
300 GB and `ingenannot_aed_time = '96h'` (four days) for future submissions.
The native command uses the requested CPU count as its worker count; speedup
depends on batch costs, memory and I/O. This does not extend an
already-running job. Extending a running Slurm time limit requires privileged
access. Native AED does not persist completed batches as resumable task outputs;
a failed task reruns its native assessment, while prepared indexed evidence can
remain cached. Changing the running process's worker count requires a new task;
resubmission restarts native assessment rather than continuing its finished batches.

`ingenannot_benchmark.nf` reads completed BRAKER and current AED outputs. It
verifies model/evidence checksums, sorts/BGZIPs/Tabix-indexes existing StringTie
GTF and miniprot GFF, runs native InGenAnnot, checks unchanged model structures
and produces matched RNA/protein comparisons. It does not download databases,
remap reads, align proteins, discover repeats, train BRAKER, refine UTRs or filter
models. Gene/transcript identities are compared exactly; changes block the report.

The root configuration records the benchmark profile:

- `ingenannot_rna_cds_only = true`
- `ingenannot_rna_penalty = 0.25`
- Protein scope is InGenAnnot's native CDS-only calculation.

These reproduce settings described in its publication, not a calibrated fish
selection policy. Our current RNA metric compares transcript exons, so differences
include coordinate scope and the native penalty as well as implementation. Native
GFF scores are rounded to four decimal places. Neither lower values nor equality
alone establishes correctness. A penalty-zero comparison can isolate that effect.

Only an environment file is added; all parameters and process resources remain
in root nextflow.config. The pinned InGenAnnot package is installed through that
Conda environment's pip section. Conda handles its dependencies. First cluster
execution verifies the environment and empty-source indexed-input behavior.

### Tiny native acceptance

After syncing, from the cluster project root:

```bash
python gcl_genome_annotation/tests/validate_ingenannot_cluster.py \
  --outdir testing/ingenannot_contracts_01
```

The test uses tiny prealigned evidence, all four source combinations, current
baseline scores, resume, profile changes and changed-evidence blocking. It stops
after native assessment; no UTR work or prediction is run. Recover failures using
the same directory and `--resume`. The original tiny native acceptance and the real benchmark passed on Crest.
The additional report-edit cache regression awaits cluster verification.

### Real unchanged-model benchmark

After tiny native acceptance, change only the Nextflow target in the existing
project sbatch to `ingenannot_benchmark.nf`. Keep the same root config, outdir,
work directory and `-resume`. It reads the completed current scoring pilot.
The initial RNA+protein pilot uses retained ETP StringTie transcripts. For later
RNA-only runs, add support for their published standalone assessment assembly
before using this development entry; absence is an explicit error, not silent
loss of RNA evidence.

Outputs are under `structural_annotation/ingenannot_benchmark/<assembly>/`:

- `benchmark_scores.tsv`: current/native RNA and protein scores, deltas and
  selected native evidence identifiers.
- `benchmark_status.json`: profile, unchanged-model checks, source provenance and
  distributions.
- `native/native.aed.gff3` and package version: original tool outputs.
- `README.md`: interpretation and links.

Inspect differences by coding/exon scope, native penalty, selected evidence and
unsupported cases. Do not introduce new annotation filters during the benchmark.

## Placement of UTR refinement

The agreed normal order is:

`BRAKER → coding baseline → UTR refinement → validate refined models → final AED/QC → model filtering → post-filter QC → functional annotation`

The existing pre-UTR baseline is retained. UTR refinement happens before final
AGAT statistics and final AED assessment. BUSCO/OMArk use proteins; if CDS/protein
identities are preserved, their scientific input is unchanged and those outputs
can be reused after verifying that invariant. Filtering can change the protein
set, requiring an appropriate post-filter comparison.

## UTR pilot

`utr_refinement.nf` now implements the next isolated development stage. Run tiny
native acceptance first:

```bash
python gcl_genome_annotation/tests/validate_utr_cluster.py \
  --outdir testing/utr_contracts_01
```

It creates two small synthetic paired BAM readsets per RNA-supported assembly
with samtools (including the multi-line native BAM-list regression), uses native InGenAnnot
0.0.20 clusterize/ranking/refinement on both strands, checks no-RNA and disabled
pass-through, extracts sequences with gffread, runs the existing AGAT module,
and verifies ranking cache reuse after a reporting-helper edit. It does not run
HISAT2, repeat discovery, BRAKER, database alignment or full-genome AED. Recover
using the same directory and `--resume`. Native cluster acceptance is pending.

After acceptance, change only the existing sbatch Nextflow target to
`gcl_genome_annotation/utr_refinement.nf`, keeping its current config, output/work
directories, resources and `-resume`. `--utr_inputs_root` defaults to `--outdir`.
The pilot reuses retained ETP StringTie transcripts or the existing standalone
RNA assembly. If neither exists for an RNA-supported assembly, it reuses the
existing StringTie assembly/merge modules on its separate accepted BAMs.

Outputs are `structural_annotation/utr/<assembly>/`: the separately validated
`<assembly>.utr.gff3`, transcript FASTA, `utr_changes.tsv`, `utr_status.json` and
README. Native ranking/refinement outputs are retained in `native/`. The existing
AGAT module also runs on the refined annotation under `qc/utr/<assembly>/agat/`,
including UTR feature types. It does not overwrite the coding baseline.

All resources and thresholds are in root nextflow.config. Ranking defaults to
16 CPUs, 64 GB and 96 hours; sequence extraction uses gffread 0.12.7. Pinned
InGenAnnot uses `clusterize -g` (not `-f`) and FOF booleans for paired and reverse
stranded BAMs. Current reverse-stranded libraries and unstranded BAMs are
supported; forward/auto strand metadata explicitly blocks the pilot until a
validated adapter is available. No BAMs are manually linked or merged.

The native tools choose RNA structures. The output adapter preserves BRAKER's
separately annotated start/stop codons and removes stop-codon bases from native
UTR features, without changing CDSs or choosing UTR boundaries itself. It records
supporting RNA IDs and strand-correct UTR lengths. The validator requires additive
exon changes only, preserved existing UTRs, gene/transcript identity, CDS/phases,
feature relationships/reference bounds and identical baseline/refined gffread
CDS/protein sequences. Original predictions remain unchanged.

The real native output contained one incompatible proposal: a minus-strand exon
under the plus-strand transcript `g10479.t2`. The adapter rejects proposals with
incompatible feature/transcript contigs or strands and restores that transcript's
original BRAKER features, including in a multi-isoform gene. Other valid proposals
remain eligible. It never flips an evidence strand to force compatibility.
`utr_status.json` records `rejected_utr_candidates`, and `utr_changes.tsv` includes
`rejected_candidate_reason`. Native outputs remain available for inspection.
Unexpected coding, identity, parent, coordinate or sequence changes still fail
validation. Restored transcripts are validated as unchanged pass-through models.
Native exon changes without a `utr_refine_evidence` tag also restore the original
transcript and are reported, rather than being accepted as supported extensions.
Candidate checks additionally reject shortening, changes to existing UTRs, and
new exonic bases inside the coding span. Native CDS coordinates, strands and
phases are checked before any restoration, so a fallback cannot conceal a coding
change. Structural and sequence validation still runs on the complete final set.

If UTR-label validation fails, inspect all disagreements before changing the
adapter or rerunning the pipeline:

```bash
python gcl_genome_annotation/tests/inspect_utr_labels.py \
  --inputs-root testing/genome_annotation \
  --assembly-id Sfl-n_primary \
  --output testing/utr_labels.json
```

This streams the original and already completed native GFFs, retains compact
coordinate structures, and reports missing/extra UTR labels after excluding
separate stop codons. Each mismatch includes transcript/evidence IDs, strand,
exon/CDS coordinates, coding-boundary equalities, and whether the candidate would
already be rejected by the adapter. It does not run ranking/refinement, read
BAMs, extract sequences or modify annotations. Inspect the otherwise-eligible
cases to determine whether the problem is label serialization or RNA structure.

The real diagnostic identified 19 missing UTR labels on 19 otherwise-eligible
transcripts. All occur where the coding boundary is a single base at the end of
an exon: InGenAnnot 0.0.20's strict inequalities omit the adjacent UTR segment.
The adapter completes only missing labels within those exact equality cases,
using native-selected exon boundaries and unchanged CDSs, excluding stop codons.
It does not extend exons, select different evidence or modify coding sequence.
Unexpected extra labels or omissions elsewhere still fail validation. Corrections
are recorded as `completed_boundary_utr_labels` in status and counted per
transcript in `utr_changes.tsv`; the unchanged native GFF remains available.

Input approval is passed to ranking as a stable content value, so editing the
validation/report helper does not invalidate the expensive native ranking when
approved scientific inputs are unchanged. This uses Nextflow's normal cache.

This pilot stops after UTR validation and AGAT. Promotion into main.nf and final
native AED rescoring are subsequent steps after cluster acceptance and review of
real extensions. BUSCO/OMArk need not be repeated for an unchanged protein set.

### AED filtering policy

AED filtering is disabled by default: `filter_max_aed = null` and
`filter_aed_source = null`. The existing downstream gene-filter module accepts
either `--filter_aed_source rna --filter_max_aed 0.7` or the corresponding
`protein` source. The threshold is a user decision; 0.7 is an example, not a
recommended calibrated cutoff. A gene is removed for AED only when every
transcript exceeds the chosen threshold; equality passes. All its isoforms are
retained when its best transcript passes, preserving the existing gene-level
filter scope. Missing evidence is an error when explicitly requested, never
silently interpreted as AED 1. Score tables retain their original columns.

The filter consumes canonical `rna_aed` / `protein_aed` columns, not benchmark
delta columns or legacy `aed`. Its caller supplies `gene_model_aed.py` as the
additional staged helper input. Other existing structural filter defaults are
unchanged. The UTR pilot does not invoke gene filtering; these choices take
effect when that downstream stage is connected after final AED assessment.

Use InGenAnnot's native short-read UTR workflow, not a custom exon-extension
algorithm: compatible merged RNA transcripts, `clusterize` against coding
models, `isoform_ranking` using existing BAMs, then `utr_refine --utr_mode rank`.
Do not take the longest assembled transcript unconditionally. Long-read evidence
can be prioritized later if supplied; current data are short-read RNA-seq.

Requirements before promotion:

1. Preserve raw BRAKER outputs and write a separate UTR-refined model set.
2. Require compatible coding structures and use native preparation and ranking.
   RNA-supported UTR/CDS overlaps can still remain after native ranking; report
   those conflicts for review rather than equating transcription with coding
   assignment or automatically deleting models.
3. Do not infer UTRs from protein-only evidence. With no RNA, skip refinement and
   retain coding models with an explicit not-applicable status.
4. Initially add UTRs to models without them, preserving existing annotated UTRs.
   Record supporting transcript IDs and 5'/3' lengths on both strands. Do not
   create arbitrary alternative coding isoforms in an additive UTR-only pilot.
5. Validate reference bounds, feature IDs/parents, CDS coordinates and phase,
   exonic containment, gene/transcript bounds and UTR orientation. Regenerate
   transcript sequences; verify CDS/protein sequence identity against the baseline.
6. Rescore refined models using native InGenAnnot and report before/after changes.
   Reuse AGAT/sequence extraction modules where compatible. Protein QC should not
   change from UTR-only work; a change is a validation failure to investigate.
7. Use raw and junction/coverage-ranked RNA support for short-read UTR assessment.
   Portcullis remains a complementary planned junction-quality assessment, not a
   hidden prerequisite or new term in native AED.
8. Keep final filtering disabled until benchmark and UTR contracts pass and real
   representative extensions are reviewed. Do not change BRAKER training inputs.

The unchanged-model benchmark, native UTR cluster contracts and real-data UTR
validation passed. Targeted BAM inspection found support for all 48 new
junctions among ten long-extension cases in at least one library (46 in both).
This supports use of the native workflow; it does not prove every assembled
transcript or its coding assignment. Final AED assessment/promotion remains
pending.

### General UTR safeguards

The wrapper preserves CDS coordinates/phases and extracted CDS/protein
sequences, IDs/parents and reference/strand consistency. Invalid optional UTR
proposals revert to the original transcript, with an explicit reason. A coding
change is a validation failure and cannot be hidden by that fallback. No-RNA
and disabled runs retain the original models.

Duplicate exon records with exactly the same coordinate set as the original
are restored to the original representation and audited separately under
`normalized_duplicate_exon_records`. They are not rejected biological
extensions. Different exon boundaries or different covered bases still require
native evidence provenance; the duplicate rule cannot accept those changes.

`utr_warnings.tsv` reports every accepted UTR addition overlapping another
gene's CDS, with strand relationship and overlapping exonic bases. Coding blocks
are merged per gene to avoid inflating counts across isoforms. Separately
annotated stop codons are excluded. Warnings never filter models or change
ranking, and the report links to them. No species-specific UTR length cutoff
is imposed. Long extensions, readthrough and overlapping coding models are
evidence conflicts to inspect, not automatic errors based on length alone.

The development utilities `review_utr_candidates.py`, `review_utr_bams.py` and
`review_utr_coding_evidence.py` inspect completed outputs without rerunning
annotation. The last selects overlapping gene pairs from the candidate report,
compares all their coding isoforms with their existing native protein AED and
selected miniprot chains, and retains protein target identifiers. The BAM
utility reports individual libraries and respects resolved strandedness. These
are targeted diagnostics, not replacements for native InGenAnnot ranking or
the planned Portcullis assessment.

## References

- [Peer-reviewed InGenAnnot study](https://pubmed.ncbi.nlm.nih.gov/40971313/)
- [Native AED command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/AED.html)
- [Native UTR workflow](https://bioger.pages-forge.inrae.fr/ingenannot/usecases/add_utrs.html)
- [UTR refinement command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/UTRRefine.html)
