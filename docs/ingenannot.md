# InGenAnnot benchmark and UTR integration

The user selected a peer-reviewed existing evaluator rather than requiring exact
MAKER equivalence. InGenAnnot 0.0.20 is the candidate. Its native AED command
provides RNA and protein scores; combined AED is omitted. The current custom
method is retained as a benchmark baseline, not silently replaced.

## First stage: unchanged-model benchmark

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
the same directory and `--resume`. Native InGenAnnot and Nextflow execution are
pending; local tests verify the validation/report adapter, not the external tool.

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

## UTR pilot requirements (next increment)

Use InGenAnnot's native short-read UTR workflow, not a custom exon-extension
algorithm: compatible merged RNA transcripts, `clusterize` against coding
models, `isoform_ranking` using existing BAMs, then `utr_refine --utr_mode rank`.
Do not take the longest assembled transcript unconditionally. Long-read evidence
can be prioritized later if supplied; current data are short-read RNA-seq.

Requirements before promotion:

1. Preserve raw BRAKER outputs and write a separate UTR-refined model set.
2. Require compatible coding structures; ambiguous, multi-gene/fusion candidates
   must be reviewed/rejected using the native preparation and ranking workflow.
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

The benchmark is implemented; the native UTR ranking/refinement integration is
the following increment after benchmark acceptance. The order and validation
requirements above do not claim UTR refinement is already implemented.

## References

- [Peer-reviewed InGenAnnot study](https://pubmed.ncbi.nlm.nih.gov/40971313/)
- [Native AED command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/AED.html)
- [Native UTR workflow](https://bioger.pages-forge.inrae.fr/ingenannot/usecases/add_utrs.html)
- [UTR refinement command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/UTRRefine.html)
