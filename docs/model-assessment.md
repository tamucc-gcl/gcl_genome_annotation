# Model assessment and repeat-derived gene flags

This stage reconnects the original pipeline's structure/sequence metrics and
TE-flagging approach to accepted post-UTR models, before model selection or
functional annotation. It is enabled in `main.nf` with `run_model_assessment=true`.
`run_te_flagging=false` disables only the repeat evidence work. No gene or
transcript is removed, renamed or altered in the accepted annotation.

## Reuse and correctness changes

- `model_metrics.py` moves sequence/exon/splice metrics from the original
  `MERGE_FUNCTIONAL_ANNOTATIONS` upstream of functional searches. The original
  GT-AG/GC-AG/AT-AC motif categories and structural severity rules are retained.
  Explicit GFF parents and unique exon/CDS blocks prevent isoform mixing and
  duplicate counting. Gene CDS/exon coverage uses union intervals across isoforms,
  not a sum of overlapping alternative transcripts. Start/stop feature presence
  and terminal protein `*` are recorded separately; their absence is not a default
  exclusion. Internal `*` detection includes the first/last residues unless the
  last symbol is the legitimate terminal marker. Functional-support and OMArk
  quality rules remain deferred until those sources are connected.
- `FLAG_TE_GENES` retains BEDTools overlap, SeqKit six-frame translation and
  DIAMOND homology, split into separately cached database/search/flag tasks.
  Repeat intervals are merged before coverage, so overlapping repeats cannot
  inflate fractions and multiple separate repeats can contribute to the 0.5
  gene-span threshold. Exon/CDS fractions provide context for intronic repeats.
- DIAMOND failures are no longer suppressed with `|| true`. A successful search
  producing an empty hit file is a valid no-hit result. Empty repeat libraries
  are reported as unavailable homology rather than manufactured negative evidence.
- Gene-level homology flags use explicit transcript parents. Transcript-specific
  hits remain visible, so one isoform's hit is not presented as a hit on every
  isoform. The derived flagged GFF preserves features, IDs, parents, coordinates,
  phases and existing notes; new attribute values are safely serialized.
- Best/worst structural severity aggregation comes from the original
  `STRUCTURAL_QC` approach. It is reported rather than used to declare a definitive
  high-confidence gene set before functional assessment.

Main uses the native classified/reclassified **pre-filter** repeat library,
consistent with the original flagging design. A supplied complete repeat library
is used when de novo discovery is bypassed. This library can contain host-gene
sequence: homology is evidence of repeat association, not definitive TE identity.
All TE flags retain models by default. Translation uses the existing standard-code
six-frame design, converts stop symbols to X and retains sequences >=30 aa.
No new identity or query-coverage removal cutoff is introduced. Query coverage is
reported from aligned query coordinates, rather than gapped alignment length.

## Outputs and semantics

Under `structural_annotation/assessment/<assembly>/`:

- `transcript_assessment.tsv`: one row per coding transcript, structural metrics
  and diagnostic `structure_flag`/reasons, separate RNA/protein AED and evidence
  IDs, junction-support columns, gene-context TE evidence and transcript-specific
  homology. Missing/disabled sources remain unavailable.
- `gene_assessment.tsv`: isoform count, best/worst structural flag, all reasons,
  protein length range and gene TE flags. No CDS sums across isoforms and no
  arbitrary combined evidence/confidence score.
- `assessment_status.json`: input identities, counts, source checksums and explicit
  zero-removal scope; `README.md` explains the results.
- `te/`: `te_gene_report.tsv`, `te_status.json`, and a derived flag-only GFF.
- `te_database/database/`: translated repeats, DIAMOND index, library identity and
  tool versions. These are project/assembly repeat data, not general shared DBs.
- `te_search/search/`: native hits, versions and protein/library/search provenance.

`FAIL` in a structural column is diagnostic severity, not a failed pipeline or
automatic model removal. Single-exon/absent-RNA junction metrics are not converted
to zero support. Native AED and junction tables are joined by transcript ID and
validated against the exact accepted model/reference identities and gene parents.
The original accepted GFF/protein/CDS sets remain the source for future selection.

The assessment, repeat overlap and reports are distinct from native homology
search. Changing only the overlap threshold reruns flags/reports; report-only
changes retain metrics and native search caches. Main passes original unchanged
coding sequence paths, avoiding UTR publication copies as new search inputs.
Changing scientific library/query/search settings correctly invalidates search.
All parameters and process resources remain in root `nextflow.config`.

## Acceptance and pilot

All 169 local Python tests pass, and the changed entry points, modules and workflows
compile with the Nextflow 23.10.1 parser. Actual SeqKit/DIAMOND/BEDTools cluster
acceptance, runtime channel wiring and real-data assessment are pending. First run on Crest:

```bash
python gcl_genome_annotation/tests/validate_assessment_cluster.py \
  --outdir testing/assessment_contracts_01
```

The test exercises actual native tools, both strands, union overlap with duplicate
repeat records, shared isoforms, four synthetic accepted evidence-mode handoffs,
valid no-hit searches, optional TE/evidence-disabled routes, main syntax, resume
and report/overlap-only cache reuse. It does not run prediction, repeat discovery,
mapping, InGenAnnot or Portcullis. Recovery uses the same directory plus `--resume`.

After cluster acceptance, use the usual main sbatch with `-resume`. An assessment-only
`model_assessment.nf` entry also reads accepted published outputs with
`--assessment_inputs_root "${outdir}"`; it requires enabled AED/junction reports.
The standalone entry finds the published pre-filter repeat library or accepts
`--assessment_repeat_library`. Missing enabled inputs are errors, never silently
disabled. Keep the same entry when testing cache behavior.

Model selection, final renaming, functional annotation and unified QC figures are
subsequent implementation increments; this stage does not activate their filters.
