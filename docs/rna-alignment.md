# RNA alignment for BRAKER3

`rna_alignment.nf` is the development entry for the next rework stage. It uses
the already accepted RNA preparation and normalized reference outputs, maps
selected reads with HISAT2 2.2.1, and produces coordinate-sorted BAMs with
samtools 1.21. It stops before BRAKER3. No repeat discovery, repeat masking,
RNA trimming, BAM merging or standalone StringTie assembly runs here.

BRAKER3 accepts multiple BAM files through `--bam=file1.bam,file2.bam` and
handles its own transcript assembly. We retain individual readset BAMs and
library metadata rather than duplicating that work. See the
[BRAKER documentation](https://github.com/Gaius-Augustus/BRAKER).

## Inputs and scope

Keep the existing sample sheet, RNA sheet, taxonomy and species arguments.
The entry defaults to the following inputs under `--outdir`:

- `rna/rna_status.json`: ready selected-read manifest from `rna_preparation.nf`.
- `filtered_assemblies/<assembly>/<assembly>.filtered.fasta` and `.contigs.json`:
  accepted normalized references from `reference_repeats.nf`.

Override their locations with `--rna_prepared_status` and `--rna_reference_dir`.
The selected FASTQ paths recorded in the RNA report must still exist. Retain
the Nextflow work directory containing processed/subset reads. Nextflow stages
these files directly, runs tasks in scratch and publishes their outputs.

This entry consumes the previously selected reads; changing trim/subset settings
requires running RNA preparation first. Changing strandedness does not require
remapping. This is a development entry, not the eventual end-to-end workflow.
When the main workflow is migrated, these stages will connect through channels.

Mapping uses the normalized unmasked reference. Repeat preparation has already
verified that the softmasked reference has exactly the same contig identities,
lengths and uppercase sequence content, so its coordinates match these BAMs.
The prediction stage will receive the accepted softmasked genome and mapped RNA.

## Alignment and validation

HISAT2 runs with `--dta`, neutral library strandedness and a reproducible seed.
It preserves distinct read groups with readset ID, library ID and sample ID.
All primary reads remain in the BAM, including unmapped reads; secondary and
supplementary records are excluded from primary-read mapping statistics.
Before mapping, the staged FASTQ byte checksums must match the accepted selected
reads, preventing an edited/replaced file from borrowing an old preparation manifest.
The BAM is sorted, indexed with CSI (including long contigs), and checked with
samtools. Full streaming validation uses pysam 0.22.1 to verify the reference
dictionary, coordinate order/bounds, read groups, layout, primary counts and
mate counts against the selected-read manifest. Zero mapped reads blocks the
readset. Mapping fraction, pairing, clipping, unique mapping and splice support
are reported; a nonzero mapping fraction alone is not a biological quality verdict.

Indexes are shared by identical reference content within the run and cached by
Nextflow. Readsets are mapped for each requested assembly. Tasks are joined by
assembly and readset identity, avoiding cross-reference associations. BAMs remain
separate even when multiple lanes belong to one library.

## Strandedness

Explicit `forward`, `reverse` or `unstranded` settings are authoritative once
BAM validation passes. `auto` uses splice-derived `XS` tags from unique primary
alignments with MAPQ >= 20. For paired data, only R1 contributes to strand
inference; both mates contribute to alignment QC and junction support.
Forward means the sequenced R1 (or single read) is in the transcript orientation;
reverse means the opposite orientation.

Default inference requires at least 100 informative reads. At least 90% in one
orientation resolves forward/reverse; 40–60% forward resolves unstranded.
Intermediate mixtures or insufficient evidence block the affected library and
appear in the report. Known library preparation settings can resolve ambiguity
explicitly. The splice-based estimate is conditional on uniquely mapped,
splice-supported reads, and may not represent all reads in an unusual library.
Conflicting lane settings within a library also block that library.

All parameters and process definitions are in `nextflow.config`:
`rna_alignment_seed`, `rna_alignment_memory`, `rna_alignment_cpus`,
`rna_alignment_threads`, `rna_no_temp_splicesite`, `rna_min_intron`, `rna_max_intron`,
`rna_strand_min_mapq`, `rna_strand_min_informative`,
`rna_strand_dominant_fraction`, and `rna_strand_unstranded_tolerance`.
Inference thresholds affect inference only; MAPQ affects QC/inference; intron
limits and alignment seed affect mapping. Unchanged tasks use native `-resume`.

Alignment defaults to 16 reserved CPUs, 14 mapper threads and 48 GB memory after
the successful mapper-only diagnostic with splice-site reuse disabled. Full
mapping with sorting at these defaults remains to be validated. CPU reservation,
mapper threads and memory are independently configurable. HISAT2 stderr also
reaches the native task error log when scratch outputs cannot be returned.
No exclusive-node flags are used. Tiny acceptance runs explicitly request four
CPUs, two mapper threads and 8 GB.

## Outputs

HISAT2 now defaults to `--no-temp-splicesite` (`rna_no_temp_splicesite = true`).
Spliced alignment and XS tags remain enabled, but splice sites discovered by
earlier reads are not reused to help later reads with short anchors. This can
affect sensitivity; mapping rate alone does not establish equivalent junction
or gene support. Set `--rna_no_temp_splicesite false` to restore reuse for a
controlled comparison. Changing this setting invalidates the alignment cache.

On 12 October 2026, full IBC1 mapping with two threads and SAM output discarded
completed with this option: 58,267,066 pairs, 76.63% overall alignment rate,
29m35s, maximum RSS 1,123,308 KiB and exit 0. The otherwise identical diagnostic
with reuse enabled grew to approximately 168 GiB before cancellation. This
strongly implicates within-run splice-site reuse in the earlier memory growth,
without establishing the underlying software defect. Full pipeline BAM/QC
validation with the new default and reduced resources remains pending.

Under `rna_alignment/`:

- `indices/<reference identity>/`: HISAT2 index, validated reference identity,
  build log and tool version.
- `readsets/<assembly>/<readset>/`: `aligned.bam`, `aligned.bam.csi`, HISAT2
  log, samtools flagstat/stats and tool versions.
- `qc/<assembly>/<readset>/`: `alignment.json` and `junctions.tsv`. Junction
  coordinates are 1-based inclusive intron coordinates; support counts alignments,
  not independent molecules.
- `strand/<assembly>/<readset>.strand.json`: supplied/inferred strand and policy.
- `README.md` and `rna_alignment_status.json`: assembly/library readiness,
  separate BAM lists, mapping statistics and reasons for blocked results.
- `multiqc/multiqc_report.html`: linked from the generated README; summarizes
  HISAT2 logs, samtools flagstat and samtools stats under assembly/readset labels.
  HISAT2's default summary format is parsed by MultiQC's Bowtie 2 / HiSAT2 module.
  The report consumes only native staged text logs, not BAMs. With no requested
  RNA, the linked page explains that no alignment was requested. Reporting uses
  MultiQC 1.25.1 and a separate helper so changes do not invalidate mapping.

The final assertion exits nonzero for incomplete evidence, after publication of
the report. Environment/scheduler/tool failures can stop execution before a report.
There is no gene prediction or transcript-model filtering in this stage.

## Cluster acceptance and real pilot

Run the tiny actual-aligner acceptance suite first, from the consuming project:

```bash
python gcl_genome_annotation/tests/validate_rna_alignment_cluster.py \
  --outdir testing/rna_alignment_contracts_01
```

The suite uses Nextflow 23.10.1 and Slurm. It creates synthetic spliced paired
and single reads and prepared manifests locally: no database downloads, repeats
or BRAKER. It checks fresh mapping, native resume, a strand-only change without
remapping, adding a distinct reference while reusing completed work, expected
blocking for insufficient strand evidence, and no-RNA routing. It verifies known
splice junctions, high mapping rates, library associations and separate BAM paths.
Cluster acceptance passed on 12 October 2026 in
`testing/rna_alignment_contracts_01` (user-supplied acceptance.json): Python
contracts, actual HISAT2 mapping recovered after an interrupted environment
build, unchanged resume, strand-only changes, adding a reference, expected
blocked exit 1 and no-RNA routing. The real alignment pilot remains pending.

An interrupted first case can be continued without rewriting its synthetic
inputs or losing the completed task cache:

```bash
python gcl_genome_annotation/tests/validate_rna_alignment_cluster.py \
  --outdir testing/rna_alignment_contracts_01 --resume
```

Stop the previous validator before resuming. Recovery completes the ready case
with native Nextflow resume, then runs the normal acceptance cases, including a
separate unchanged-resume check. Environment creation failures are infrastructure
failures and do not imply that alignment failed.

After acceptance, change only the Nextflow script in the personal sbatch command
to `gcl_genome_annotation/rna_alignment.nf`. Keep the existing input/taxonomy
arguments, `--outdir /work/birdlab/GCL/prj_seemann-flounder_genome/testing/genome_annotation`,
work directory and `-resume`. This uses the completed real RNA/reference outputs,
adds the alignment results, and stops before BRAKER3. No need to cancel it to
avoid rerunning repeats: repeat processes are absent from this entry.

Alignment settings follow the [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/).
