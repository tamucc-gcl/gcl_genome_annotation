# gcl_genome_annotation
Pipeline for Genome Annotation

![](assets/genome_annotation_pipeline.png)
## Refactor foundation (cluster validation pending)

The new [input and taxonomy preflight](docs/foundation-preflight.md) supports two-sheet validation, global defaults, per-assembly overrides and taxonomy-derived settings. It is a separate entry point; the legacy annotation workflow has not yet been migrated. See the handoff for cluster commands and current limitations.

The [rework addendum](docs/rework-decisions.md) records
end-to-end execution with native resume and existing model-filter thresholds as
the default evaluation baseline.
