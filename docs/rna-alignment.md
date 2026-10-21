# RNA alignment and strand QC

The integrated workflow maps accepted selected reads to normalized unmasked
references using HISAT2 2.2.1 and samtools 1.21. BAMs remain separate by assembly
and readset, including lanes of the same library. There is no pipeline BAM merge.
`rna_alignment.nf` is a diagnostic entry using completed preparation/reference
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
resolves unstranded. Intermediate/insufficient data block inference, with known
protocol overrides available. Conflicting lane strands within a library block it.
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
with the softmasked prediction reference. Strand-only changes need not remap;
changed reads/reference/mapping parameters correctly invalidate dependent tasks.

```bash
python gcl_genome_annotation/tests/validate_rna_alignment_cluster.py \
  --outdir testing/rna_alignment_contracts_01
```

Acceptance exercised actual spliced paired/single mapping, separate BAMs, known
junctions, recovery/resume, strand-only changes, added references, blocked
inference and no RNA. The full paired-end real-data mapping/QC pilot completed.
See [running](running.md) and [BRAKER](braker-prediction.md).
