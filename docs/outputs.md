# Output directory guide

Start with `<outdir>/annotation_report.md`, then `release/README.md` and each
assembly's release README. Stage-native/audit outputs remain available separately.

| Directory | Current contents |
|---|---|
| `preflight/` | Resolved inputs, taxonomy, origins and readiness |
| `protein_sources/`, `protein_evidence/` | Source provenance and per-assembly combined prediction evidence |
| `rna/` | Raw/processed/selected read manifests, FastQC/fastp, separate raw and filtered MultiQC |
| `filtered_assemblies/<assembly>/` | Normalized retained reference, contig identities/header map and exclusions |
| `reference_repeats/`, `repeats/` | Repeat preparation status, libraries, masking outputs, union coverage, logs/audits and native landscape |
| `rna_alignment/` | HISAT2 indices, separate BAM/CSI pairs, mapping/strand QC and MultiQC |
| `structural_annotation/braker_inputs/`, `braker_runtime/` | Prediction handoff and software/runtime records |
| `structural_annotation/braker/<assembly>/` | Original BRAKER predictions, coding sequences, logs and trained parameters |
| `structural_annotation/utr/<assembly>/` | Native UTR candidates, accepted changes/rejections, warnings and unfiltered extracted transcripts |
| `structural_annotation/final/<assembly>/` | Accepted post-UTR, pre-selection models and coding FASTAs; despite its name, this is not the exported final selected release |
| `structural_annotation/aed/<assembly>/` | Native RNA/protein AED and evidence provenance |
| `qc/junctions/<assembly>/` | Portcullis native/readset results and model junction support |
| `structural_annotation/assessment/<assembly>/` | Complete structural, sequence, AED, junction and TE evidence assessment |
| `structural_annotation/selection/<assembly>/` | Retained/excluded sets and every selection decision |
| `structural_annotation/identified/<assembly>/` | Selected final IDs, complete maps and re-keyed assessment |
| `qc/final/`, `qc/selected/` | Unfiltered post-UTR versus selected AGAT/BUSCO/OMArk; representative FASTAs/tables |
| `qc/structural/` | Markdown before/after structural reports, complete tables, ten PNG/SVG figure groups |
| `functional_annotation/<assembly>/` | Native functional results, eligibility, merged tables and decorated coding GFF3 |
| `genome_tracks/<assembly>/` | Default selected tracks, raw tracks, filter policies/decisions and status |
| `ncRNA/`, `telomeres/`, `centromere/` | Native track outputs; these retain raw search meaning |
| `release/<assembly>/` | Copied, checksummed current deliverables; see below |
| `reports/` | Unified report data snapshots, metrics, plots and report status |
| `pipeline/` | Native Nextflow trace, execution report, timeline and DAG |

Standalone QC/benchmark entries can additionally publish `qc/baseline/`,
`structural_annotation/ingenannot_benchmark/` and custom evidence-score outputs.
Those directories are not production stages merely because an earlier run left
them behind. The current [entry-point inventory](testing.md) distinguishes them.

## Release payload

Each release includes combined and coding-only GFF3, unmasked and softmasked
references, selected protein/CDS/spliced-transcript FASTAs, an enabled selected
tRNA FASTA, complete ID maps, assessment/functional tables, separate tracks,
provenance, `release_status.json` and `checksums.sha256`. The transcript FASTA
includes accepted UTRs; tRNAs are intron-spliced. Both use final selected IDs.
Enabled no-hit tRNAs yield an empty FASTA; disabled tRNAs yield none.

```bash
cd /path/to/outdir/release/ASSEMBLY_ID
sha256sum -c checksums.sha256
```

The current manifest defines deliverables. Older published files may persist
after settings change; presence alone does not make them current. Keep the output
hierarchy together for report links. Moving just the release preserves its
payload, but links back to detailed stage reports require the surrounding tree.
Intermediate linked files still depend on original inputs/work; copied releases
do not. See [packaging safeguards](annotation-release.md).
