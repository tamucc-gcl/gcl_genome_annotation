# Documentation

These guides describe the connected `main.nf` workflow. The root
[configuration](../nextflow.config) is authoritative for defaults and process
resources. Historical plans are clearly separated from current instructions.

## Operation

- [Running, databases, resources and resume](running.md)
- [Input sheets, overrides and taxonomy](foundation-preflight.md)
- [Output directory and deliverables](outputs.md)
- [Testing and entry points](testing.md)
- [Validation record and limits](validation.md)

## Stages

| Stage | Guide |
|---|---|
| Protein sources | [Protein evidence](protein-evidence.md) |
| RNA preparation and QC | [RNA preparation](rna-preparation.md) |
| References and repeats | [Reference/repeat preparation](reference-repeats.md) |
| RNA mapping and strand | [RNA alignment](rna-alignment.md) |
| Prediction | [BRAKER3](braker-prediction.md) |
| UTR refinement | [InGenAnnot](ingenannot.md) |
| Evidence support | [Native RNA/protein AED](evidence-scoring.md) |
| Splice support | [Portcullis](junction-support.md) |
| Unfiltered QC | [AGAT/BUSCO/OMArk baseline](annotation-baseline.md) |
| Assessment and selection | [Assessment](model-assessment.md), [selection](model-selection.md) |
| Final IDs and before/after QC | [Final structural annotation](final-structural-qc.md) |
| Functional searches | [Functional annotation](functional-annotation.md) |
| Independent features | [Genome tracks and their filters](genome-tracks.md) |
| Delivery and report | [Release](annotation-release.md), [unified report](unified-report.md) |

## Maintenance and history

[Cleanup decisions and file inventory](cleanup.md) record completed source consolidation
and remaining archival candidates. [The original/current audit](rework-audit.md) records the comparison
at its stated commit and subsequent sequence-export follow-up.

The following are historical decision/design records, not instructions to run a
second workflow: [rework decisions](rework-decisions.md), [old continuation plan](evidence-scoring-plan.md),
[MAKER compatibility proposal](aed-requirements.md), and [custom-score evaluation](evidence-scoring-evaluation.md).
The implemented method is native InGenAnnot with separate RNA/protein AED.
