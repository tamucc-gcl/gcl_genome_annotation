# Unfiltered and selected scientific QC

`ANNOTATION_BASELINE` reuses AGAT, BUSCO and OMArk to assess accepted annotations.
In `main.nf`, unfiltered QC follows UTR validation under `qc/final/`; selected
final-ID QC follows selection under `qc/selected/`. Original BRAKER results remain
preserved. The focused `annotation_qc.nf` entry instead assesses raw accepted
BRAKER output under `qc/baseline/` and stops without prediction or filtering.

## Populations and tools

AGAT assesses all GFF3 isoforms, including UTRs, and retains native statistics,
YAML and distribution outputs. BUSCO 5.3.2 runs in **protein mode** with the
resolved lineage. BUSCO and OMArk receive one longest protein per coding gene,
not the genome FASTA or all isoforms. Terminal `*` is excluded from length
ranking, not silently removed from the scientific input sequence.

Ties use transcript IDs; selected QC preserves the original-ID tie-break through
renaming. Exact representative FASTAs/tables remain available. Selection of a
QC representative does not remove genes or replace the complete proteome.
Per-assembly `run_busco`/`run_omark` switches are honored. Disabled analyses are
`not_requested`, not zero completeness. Native unclassified/short-protein
behavior remains explicit; unassessed isoforms do not inherit OMArk classes.

Shared stores retain `busco_database/lineages/<lineage>/` and
`omark_database/<oma_db_name>.h5` plus `taxa.sqlite`. Existing contents are
reused natively; declared settings are not independent database-release audits.
Full native BUSCO/OMArk directories, summaries and versions are preserved.

## Controls and outputs

`run_baseline_qc=true` enables unfiltered main QC. Selected QC is part of
`run_finalization`; per-tool BUSCO/OMArk switches control both populations.
Each QC root contains its README/status and per-assembly `selection/`, `agat/`,
`busco/`, `omark/` outputs as applicable. A QC input named `<assembly>.fa` here
is the representative **protein** FASTA, not the genome. Ready indicates
requested tool completion and integrity, not a completeness threshold.

```bash
python gcl_genome_annotation/tests/validate_annotation_baseline_cluster.py \
  --outdir testing/annotation_baseline_contracts_01
```

This focused acceptance uses actual AGAT, representative selection and resume/
integrity cases; it does not run BUSCO/OMArk's full databases. Their real-data
unfiltered and selected runs were reviewed separately. The initial raw BRAKER
pilot and later integrated model counts differ; use current manifests for a
current comparison, not historical pilot totals. See
[final structural QC](final-structural-qc.md) and [validation](validation.md).
