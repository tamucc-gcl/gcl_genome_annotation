# Unified Markdown annotation report

The approved format is **one complete Markdown report**, with an overview and
assembly-specific subsections. Key figures are inline; detailed figures and tables
are linked. Content is restricted to results and methods: no tuning advice,
interpretation, composite confidence score or pipeline-written HTML.

Open `annotation_report.md` at the output root. Supporting outputs are:

- `reports/annotation_metrics.tsv`: explicit per-assembly metrics and denominators.
- `reports/annotation_report_plots/`: new overview/detail PNG and SVG figures.
- `reports/annotation_report_data/`: source stage-manifest snapshots.
- `reports/annotation_report_status.json`: inputs, output checksums and scope.

The report retains the original Markdown table/number-formatting helpers and
section subjects in `generate_annotation_report.py`. The production entry now
uses `unified_report.py` and the accepted keyed manifests. Legacy filename/source
assumptions (`BRAKER3`, `g123.t1`), obsolete output links, combined/legacy AED and
threshold-based tuning suggestions are removed from production reporting.

The old R plotting script remains available as a reference. Current structural
figures already provide accepted before/after structure, completeness, AED and
Portcullis views. These are reused, not recomputed. The report carries forward
functional coverage, functional-vocabulary/product frequency and telomere-position
plots using current release IDs and the existing Matplotlib environment. Vocabulary
plots are frequencies, not enrichment analyses; terms are deduplicated within each
selected transcript. Product frequencies count coding genes.

## Per-assembly layout

1. Release links and reference/evidence inputs.
2. BRAKER, optional UTR changes and selection counts/settings.
3. Before/after completeness and structural distributions.
4. Separate RNA/protein AED and coding/UTR-associated junction support.
5. Repeat coverage, raw/selected genome tracks and filter decisions.
6. Functional products, eligibility and tool coverage.
7. Methods, resolved settings, stage statuses and source reports.

The inline structural figures are retention, completeness, separate AED and
splice-support panels. New inline figures cover track filtering, selected telomere
positions and functional coverage when applicable. Other structural and detailed
functional figures are linked. Native MultiQC and other tool HTML remain linked.

## Population and provenance checks

Release counts are authoritative. Repeat records, ncRNA matches, coding genes and
transcripts remain distinct populations. RepeatMasker RNA-labelled repeats are
not added to ncRNA gene counts. Telomere window decisions and merged regions
remain distinct. BUSCO/OMArk use representative proteins; junction occurrences
across isoforms differ from distinct genomic junctions. Native RNA/protein AED
remain separate and missing evidence is not replaced with a numeric score.

Configured versions and selectors (including `CURRENT`) are labelled as settings,
not proof of the version in a pre-existing shared store. Native runtime records
and source provenance remain linked. The report verifies consumed release/QC
checksums, reference/model identities and complete per-assembly stage handoffs.
Failed or missing required inputs block reporting; they cannot become zeros or
an all-clear statement. No additional scientific filtering is performed.

`main.nf` invokes reporting after validated release packaging and the required
source reports. Report-helper edits rerun only reporting. Existing publication
paths and release payloads are unchanged. Keep the output hierarchy together to
resolve relative links; plotting/report outputs are published as copies.

## Acceptance

After syncing, run on Crest with the usual Nextflow/Conda modules:

```bash
python gcl_genome_annotation/tests/validate_report_cluster.py \
  --outdir testing/report_contracts_01
```

This runs actual Nextflow/Matplotlib reporting for two accepted synthetic releases,
including mixed evidence modes, enabled empty/disabled tracks, current counts,
source checksums, PNG/SVG output, unchanged resume and report-only reruns. It
does not run scientific searches. Use `--resume` with the same directory after
a failure. Synthetic source structural plots are placeholders; native structural
plot correctness was validated by the earlier structural-QC stage tests.

After acceptance, resume the usual `main.nf` run. Review `annotation_report.md`
and `reports/annotation_report_status.json`. Native cluster reporting acceptance
and the full-genome report review remain pending.
