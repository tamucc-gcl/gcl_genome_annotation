# RNA alignment and strand QC

The integrated workflow maps accepted selected reads to normalized unmasked
references using HISAT2 2.2.1 and samtools 1.21. BAMs remain separate by assembly
and readset, including lanes of the same library. There is no pipeline BAM merge.
`validation/entrypoints/rna_alignment.nf` is a diagnostic entry using completed preparation/reference
outputs; it stops before BRAKER. Main uses native task channels throughout.

HISAT2 uses `--dta`, a fixed seed, strand-neutral mapping and explicit read-group
IDs/library/sample. Default intron limits are 20..500,000 bp. Before mapping,
FASTQ byte checksums must match the selected-read manifest. The BAM is coordinate
sorted, CSI-indexed and checked by samtools; streaming pysam validation checks
reference dictionary, bounds/order, read groups, layout and primary/mate counts.
Zero mapping blocks a readset. Primary statistics exclude secondary/supplementary
records; alignment percentages alone are not a biological accuracy verdict.

## Strand resolution

Explicit forward/reverse/unstranded settings are retained after validation.
`auto` uses unique primary, splice-derived XS-tag observations with MAPQ >=20;
paired libraries use R1 for inference. At least 100 informative reads are needed.
At least 90% in one orientation resolves forward/reverse; 40..60% forward
resolves unstranded. These thresholds have not changed. Intermediate proportions
remain `ambiguous`; fewer informative observations remain `insufficient`. Both
retain `inferred_strandedness=unresolved`.

`rna_unresolved_strand_policy='unstranded'` is the default in main and the
alignment diagnostic entry. Otherwise valid unresolved BAMs proceed with no
library-strand assumption and visible warnings. `--rna_unresolved_strand_policy error`
requires confident automatic inference and blocks unresolved libraries with an
actionable reason. Invalid policy values are configuration errors. Neither policy
relaxes BAM integrity, zero-mapping, identity, layout or required splice checks.

For automatic libraries, disagreement or uncertainty in any lane selects
unstranded processing for **all lanes of that physical library**. Every lane's
inference and measurements remain intact. Strict policy blocks disagreement.
Explicit settings retain their source and observed inference evidence; inconsistent
explicit settings or mixed requested settings in a library remain configuration
errors at preflight. Fallback never overrides explicit instructions.
Forward means R1/single-read orientation matches the transcript; reverse is the
opposite. This splice-conditional estimate need not represent every RNA molecule.

Default reservation is 16 CPUs/48 GB, using 14 mapper threads. Within-run temporary
splice-site reuse is disabled (`rna_no_temp_splicesite=true`) following the
observed memory growth in a real diagnostic. Spliced mapping and XS tags remain;
sensitivity to short anchors can change. Changing this option invalidates mapping.
The corrected full pipeline BAM/QC run and MultiQC report were accepted. Resource
requests are configurable and are not guaranteed runtime memory caps.

## Outputs and reuse

`rna_alignment/` contains reference-content indices; separate
`readsets/<assembly>/<readset>/aligned.bam` and `.csi`; mapping/strand QC,
native HISAT2/samtools logs, `rna_alignment_status.json`, README and
`multiqc/multiqc_report.html`. MultiQC consumes text logs, not staged BAM copies.
No-RNA results remain explicit. Preserve BAM inputs and backing work when linked.

The focused entry reads `rna/rna_status.json` and `filtered_assemblies/` under
outdir, overridable with `rna_prepared_status`/`rna_reference_dir`. It does not
rerun preparation when trim/subset settings change. Reference identity is shared
with the softmasked prediction reference. Policy and explicit strand values are excluded from mapping job specifications.
The independent `rna_mapping.py` helper isolates mapping/reference checks from
strand inference/report edits. Within this revision, policy-only changes should
reuse indices and BAMs with `-resume` in the same entry/work directory, while
strand analysis, reports and affected downstream tasks are invalidated. The helper
split changes task inputs once on upgrade: old cached index/mapping tasks may
rerun. Do not assume all downstream results or caches across entrypoints can be
reused. Changed reads/reference/mapping parameters invalidate dependent tasks.

```bash
python gcl_genome_annotation/tests/validate_rna_alignment_cluster.py \
  --outdir testing/rna_alignment_contracts_01
```

Acceptance exercised actual spliced paired/single mapping, separate BAMs, known
junctions, recovery/resume, strand-only changes, added references, strict blocked
inference and no RNA in the earlier implementation. New default-fallback/strict
policy resume cases are added to the controller and require cluster execution. The full paired-end real-data mapping/QC pilot completed.
See [running](running.md) and [BRAKER](braker-prediction.md).

## Machine-readable contract (schema version 2)

`rna_alignment_status.json` and per-readset strand manifests distinguish:

| Field | Meaning |
|---|---|
| Readset `status`, strand `alignment_status` | Hard alignment QC, independent of inference |
| `requested_strandedness` (`requested` alias) | auto / forward / reverse / unstranded |
| `inferred_strandedness` | forward / reverse / unstranded / unresolved |
| `inference_outcome`, `inference_reason` | confident / ambiguous / insufficient / not_evaluated (BAM failure), and reason |
| `effective_strandedness` | Operational setting, null when no setting is accepted |
| `effective_setting_source` | explicit / confident_inference / conservative_fallback |
| `informative`, `strand_forward`, `strand_reverse`, fractions | Unmodified observations, including with explicit settings |
| `method`, `policy` | splice_XS and thresholds/requested unresolved policy |
| `warning`, top-level `warnings`, `has_warnings` | Uncertainty and selected fallback with assembly/library/readset identities and measurements |

The compatibility field `strandedness` means **effective processing** in v2;
it must never be used to describe inferred preparation. Library rows contain the
agreed processing setting; the aggregate report applies it to every member
readset while retaining each lane's inference. Per-readset strand files describe
the lane decision before aggregation; downstream consumers use the aggregate.
Ready with warnings means operationally usable, not warning-free or biologically
validated. README and the unified report retain this distinction.

`spliced_xs_reads` counts mapped primary spliced alignments with a usable `XS:+/-`
tag across both mates. BRAKER requires at least one such record per requested BAM
in addition to existing usable-splice checks, without requiring tags on every
alignment. Legacy QC without this measurement must be rerun with current helpers;
do not edit status JSON to upgrade it. `samtools`/pysam integrity/index checks remain.

StringTie uses `--fr` for forward, `--rf` for reverse and neither for unstranded
processing, including fallback. HISAT2 remains strand-neutral (`--dta`); XS tags
are never rewritten. Portcullis uses native `unstranded`, independently of FR/RF
mate orientation. InGenAnnot uses its supported unstranded ranking mode; known
forward inputs still require disabling UTR refinement or a future validated
adapter. BRAKER receives validated BAMs in ET/ETP, with neither `--stranded` nor
experimental UTR prediction enabled. Retained ETP StringTie transcripts come
from that strand-neutral route; otherwise AED/UTR transcripts are reassembled
using effective settings. See [BRAKER](braker-prediction.md),
[Portcullis](junction-support.md), and [InGenAnnot](ingenannot.md).
