# RNA acquisition and preparation (synthetic cluster acceptance passed)

The real protein pilot and its stored-source rerun completed on 9 October 2026.
The next increment is `rna_preparation.nf`: input preflight, acquisition, FASTQ
integrity checks, raw QC, trimming, processed QC, optional subset and reporting.
It stops automatically before reference alignment, repeat work or BRAKER, and does
not rerun protein preparation. These temporary development entry points will be
connected in the normal uninterrupted workflow as the remaining stages are built.
All configuration is in the root `nextflow.config`.

## Inputs and identity

Keep using the existing assembly and RNA CSV files and their readset associations.
Only readsets requested by at least one assembly are prepared; unused sheet rows
are still checked by configuration preflight. Each requested readset is processed
once, even if multiple assemblies use it. Readsets and physical libraries stay
distinct; lanes are not concatenated, and later alignment can retain their groups.
Local FASTQ and gzip FASTQ are supported for both paired and single-end inputs.

Run accessions use pinned SRA Toolkit 3.1.1: `prefetch` followed by `fasterq-dump
--split-files`. Preserve accession requests, acquisition commands/logs, downloaded
SRA files, compressed raw FASTQ and tool version. Single-end must produce exactly
one expected read file; paired-end must produce exactly the requested two mates.
Full-content validation then checks the same route used for local inputs.
Uncompressed conversion output and conversion scratch are temporary; they are
removed after successful raw validation. Allow enough disk for SRA conversion:
[NCBI documents its temporary space requirements](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump).
Real network accession acquisition has not yet been tested; unit fixtures mock it.

FASTQ validation scans every record, including decompression integrity, four-line
structure, nucleotide alphabet, sequence/quality length and printable quality
characters. Paired files must have equal counts and matching ordered IDs;
`/1`, `/2` and supplied CASAVA mate fields must agree with mate order. Empty input
and empty selected output fail. Phred+33 is the expected input encoding; printable
range checks cannot independently distinguish all Phred+64 data. Wrapped FASTQ
records are outside the initial input format.

## Processing defaults

FastQC 0.12.1 runs on full raw reads and full cleaned reads. Names include the
readset, raw/processed stage and mate so results remain distinct in MultiQC 1.25.1.
Native HTML/ZIP reports, FastQC commands/versions, fastp JSON/HTML/commands/logs,
read retention and input/output SHA256 identities are preserved. QC warning/fail
flags are reported; they do not automatically discard a readset.

[fastp 0.24.0](https://github.com/OpenGene/fastp) performs adapter removal and
conservative 3-prime tail trimming when enabled:

| Setting | Default |
|---|---|
| `trim_adapters` | true; paired-end adapter auto-detection enabled |
| `trim_quality` | true; `cut_tail`, four-base window, Q20 mean |
| `rna_trim_quality` | 20 |
| `rna_trim_window` | 4 |
| `rna_min_read_length` | 30 bases after trimming |
| `rna_max_n_bases` | 5 uncalled bases per read |
| `subset_fraction` | 1.0, all cleaned reads |
| `subset_seed` | 42 |

RNA-sheet boolean overrides take precedence over global defaults. With both
trimming flags false, consumers reuse the existing reads and no fastp filtering
runs. With either enabled, the configured minimum-length and maximum-N
filter apply. General percentage-of-low-quality filtering and automatic polyG
trimming are explicitly disabled. No deduplication, merging, orphan-read rescue,
polyA trimming or low-complexity filtering is added. Report these initial defaults
and evaluate them on the real RNA before changing thresholds.

The optional subset follows cleaned-read QC and uses a deterministic SHA256
selection on seed and read ID. Both mates are selected together. The fraction is
an expected fraction rather than an exact count; duplicate read IDs receive the
same decision. Full cleaned reads are retained regardless of subset fraction.
Fraction 1 passes the existing read paths onward with metadata only. A subset producing no
reads blocks its consumers rather than silently supplying empty alignment input.

Strandedness is preserved as supplied (`unstranded`, `forward`, `reverse`) or
`auto`. Auto remains **pending reference-specific alignment**. This stage does
not infer strand or silently turn auto into unstranded, and its preparation-ready
status is not permission to bypass strand inference in the future prediction stage.

## Outputs and failures

```text
<outdir>/rna/
    README.md
    rna_status.json
    raw/<readset_id>/rna.json; downloaded reads/source logs for accessions only
    processed/<readset_id>/rna.json; newly cleaned reads/fastp reports when trimming
    selected/<readset_id>/rna.json; newly selected reads when fraction < 1
    qc/raw/<readset_id>_raw/ native FastQC reports
    qc/processed/<readset_id>_processed/ native FastQC reports
    multiqc/raw/multiqc_report.html and data/logs (before filtering)
    multiqc/filtered/multiqc_report.html and data/logs (after filtering + fastp)
```

Newly produced read bundles are published by Nextflow as links to retained work, avoiding
copies for every consuming assembly. Preserve that work/cache and original local
FASTQ inputs. RNA data is project evidence, not a general database; it uses native
`-resume`, rather than the general-database shared `storeDir` policy.
Two MultiQC reports separate raw FastQC from processed FastQC; the processed report
also includes fastp statistics when trimming ran. Explicitly requested modules use
`--require-logs`, so missing expected module results fail report generation rather
than silently disappearing. No-RNA runs emit a clearly labelled no-input page for
each stage. These reports describe full raw/cleaned reads before optional subsetting.
Local FASTQs remain direct path inputs. Helpers do not create pass-through links
or read copies. FastQC input names are assigned through Nextflow staging, and the
report stages small manifests/QC outputs rather than read bundles. Read paths in
`rna_status.json` identify the actual selected inputs; preserve their source/work
locations. Processing tasks inherit scratch and staging settings from the root
configuration. See [Nextflow input staging](https://docs.seqera.io/nextflow/process).
Expected acquisition/content/trimming/subset failures become explicit readset
failures. Dependent assemblies are blocked, independent readsets can finish, and
the final assertion exits non-success after reporting. QC-tool or infrastructure
failures remain task failures and can stop the workflow before final reporting.

Raw acquisition descriptors exclude taxonomy, assembly lists, trimming, subset and
strand settings. Trimming descriptors exclude subset and strand settings. Subset
descriptors depend on cleaned reads plus fraction/seed. Changing subset settings
should reuse full-read acquisition, raw QC, trimming and processed QC. Report
logic is separate so report edits need not reprocess reads.

## Cluster acceptance and real pilot

Implementation is source-reviewed, with Python syntax checks passed. Tests and
Nextflow execution are reserved for the cluster. After syncing, run:

```bash
python gcl_genome_annotation/tests/validate_rna_cluster.py \
  --outdir testing/rna_contracts_01
cat testing/rna_contracts_01/acceptance.json
```

The runner executes the unit suite plus tiny actual-tool paired/single fixtures,
shared RNA consumers, an unchanged resume, a subset-only settings change, a shared
bad-mate failure with an independent successful assembly, and no-RNA routing.
It downloads no real biological data; worker environments are created as needed.
The real available pilot is paired-end, so successful single-end fixtures do not
constitute biological single-end validation.

For the real pilot, change only the entry script in the existing sbatch command:

```bash
run "${PROJECT_ROOT}/gcl_genome_annotation/rna_preparation.nf"
```

Keep the single root config, sheets, project/database arguments, testing output
and log directories, and `-resume`. Both real pairs already associated in the RNA
sheet can run; each will be prepared once. Review `rna/README.md`, the MultiQC
report and retention numbers. This run stops after RNA preparation. Next comes
reference/repeat integrity, followed by reference-specific alignment and strand
resolution before prediction integration.

First cluster attempt on 9 October 2026: the unit suite passed, then Nextflow stopped while creating the acquisition Conda environment before publishing the RNA report. Local inputs now use Python alone; only run-accession inputs install SRA Toolkit. RNA tool environments use version pins with conda-forge ahead of bioconda, and explicit channel overrides exclude inherited/default channels. See [Bioconda setup](https://bioconda.github.io/). The corrected cluster solve and actual-tool acceptance remain pending; use testing/rna_contracts_02 for the next attempt.

Cluster acceptance passed on 9 October 2026 in testing/rna_contracts_02 (user-supplied acceptance.json): unit tests; fresh paired/single/shared-readset preparation with actual QC/trimming tools; unchanged resume; subset_fraction=0.5 changing only selection; shared failed-readset isolation (expected exit 1); and no-RNA routing. The corrected Conda environments solved and ran successfully. Real paired-end RNA preparation and real network accession acquisition remain pending. Automatic strandedness remains deferred to reference-specific alignment.

Subsequent cleanup removes helper-created local/passthrough links and stages actual
FASTQs as Nextflow path inputs. Processing tasks inherit configured scratch/staging;
reporting stages metadata rather than FASTQs. The first cleanup acceptance attempt
stopped at module parsing; dynamic input path declarations corrected that error.
The corrected cleanup passed cluster acceptance on 9 October 2026 in
`testing/rna_contracts_04` (user-supplied acceptance.json): unit tests, fresh paired/
single/shared processing, unchanged resume, subset-only reruns, disabled-trim and
all-read passthrough fresh/resume, expected blocked failure (exit 1), and no-RNA
routing. This validates RNA native staging and configured scratch operation;
the separate changed protein handoffs still require protein acceptance. Real RNA
preparation, real accession downloads, and reference-specific strand inference
remain pending.

Accepted cluster command:

```bash
python gcl_genome_annotation/tests/validate_rna_cluster.py \
  --outdir testing/rna_contracts_04
```

The real local pilot subsequently reported ready for both readsets: 58,267,066
selected pairs for IBC1_S123_L006 and 58,769,470 for IBC2_S124_L006. The uploaded
combined MultiQC HTML contains all eight raw/processed FastQC samples and fastp
results; the user confirmed it displays correctly. Report generation now produces
separate raw/filtered pages. Python syntax checks passed; the split reports and
updated report-content acceptance checks await cluster validation. Resume the
existing real run to regenerate reports from cached processing, or run the cluster
acceptance script with a new empty output directory such as testing/rna_contracts_05.
