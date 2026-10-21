# gcl_genome_annotation

A Nextflow DSL2 pipeline for eukaryotic nuclear genome annotation. The main
workflow prepares references and evidence, predicts coding genes with BRAKER3,
assesses and selects models, annotates their functions and independent genome
tracks, and packages a validated release with a unified Markdown report.

## Start here

- [Run, configure and resume the pipeline](docs/running.md)
- [Input sheets and taxonomy](docs/foundation-preflight.md)
- [Outputs and final deliverables](docs/outputs.md)
- [Tests and validation scope](docs/testing.md)
- [Documentation index](docs/index.md)

Run from the consuming project directory, with this repository checked out as
`gcl_genome_annotation/`. The root [nextflow.config](nextflow.config) holds all
parameter defaults, tool environments, resources and profiles. Review its
site-specific paths and scheduler settings before using another cluster.

```bash
nextflow -C gcl_genome_annotation/nextflow.config \
  run gcl_genome_annotation/main.nf \
  -profile slurm -resume -w work \
  --sample_sheet data/annotation_assemblies.csv \
  --rna_sheet data/annotation_rna.csv \
  --taxid YOUR_TAXID --taxdump_dir /path/to/ncbi_taxonomy \
  --dfam_database /path/to/databases/dfam \
  --diamond_database /path/to/databases/diamond \
  --orthodb_database /path/to/databases/orthoDB \
  --rfam_database /path/to/databases/rfam \
  --busco_database /path/to/databases/busco \
  --omark_database /path/to/databases/omark \
  --interproscan_database /path/to/databases/interproscan \
  --eggnog_database /path/to/databases/eggnog \
  --kofamscan_database /path/to/databases/kofamscan \
  --outdir results/genome_annotation
```

Replace the example paths and taxid. Omit `--rna_sheet` when no assembly uses
RNA. Disable unused analyses/evidence explicitly; omitting a database path does
not disable its analysis. See the running guide for existing database layouts,
mode selection, dependencies and resource controls.

The production entry is `main.nf`. Focused stage launchers are kept under
`validation/entrypoints/` for diagnostics and regression testing; see the
[testing guide](docs/testing.md). Development scratch material is not required
for running the pipeline.

## Workflow

| Order / branch | Work |
|---|---|
| Inputs | Validate assembly/RNA sheets; resolve taxonomy and per-assembly settings |
| Preparation | Normalize/filter references; acquire/combine proteins; validate, trim and QC RNA |
| Repeats | Dfam, RepeatModeler, curation/classification/Terrier, protein-based library filtering, parallel RepeatMasker |
| Mapping and prediction | HISAT2 with separate readset BAMs; validated BRAKER3 evidence routing |
| Postprocessing | Optional native InGenAnnot UTR refinement, unfiltered QC, native RNA/protein AED and Portcullis junction assessment |
| Selection and final IDs | Structural/evidence rules, retained/excluded audits, deterministic IDs and complete maps |
| Selected QC and function | AGAT/BUSCO/OMArk before/after comparison; InterProScan, eggNOG-mapper, DIAMOND, KofamScan and DeepSig |
| Independent tracks | tRNAscan-SE, Infernal/Rfam, tidk and centroAnno; raw/selected QC audits |
| Delivery | Combined/separate GFF3s, references, protein/CDS/transcript/tRNA FASTAs, tables, maps and checksums; unified Markdown report |

Independent branches run concurrently when their inputs are available. This is
one normal end-to-end invocation; standalone entries exist for focused tests and
diagnostics, not required manual handoffs.

## Evidence and selection

| Requested evidence | BRAKER route |
|---|---|
| RNA + proteins | ETP |
| RNA only | ET |
| Proteins only | EP |
| Neither | ES self-training |

RNA is associated by explicit readset ID. Protein evidence can combine OrthoDB,
NCBI assembly accessions and local FASTAs. The planner blocks missing requested
evidence instead of silently downgrading modes. `braker_ab_initio=false` controls
extra no-hint AUGUSTUS output, not whether genome-only ES prediction can run.

The default structural selection removes single-exon transcripts with CDS
<300 bp, transcripts in genes spanning <200 bp, and spliced transcripts whose
introns are all <10 bp. A gene survives if any isoform survives. AED and junction
filtering are off by default; TE association is flag-only. All decisions and
excluded models remain available. [Selection details](docs/model-selection.md).

## Results

Open `results/genome_annotation/annotation_report.md` for the complete report,
then `release/README.md` and `release/<assembly>/README.md` for deliverables.
The report uses assembly subsections, key inline figures and links to detailed
tables/native results. Pipeline-written reports are Markdown; native tool HTML
reports remain. The current release manifest is authoritative if older published
files remain after settings change.

Nextflow manages staging, scratch, transfers and task caching. Keep the same
launch directory, `.nextflow/` cache, work directory and inputs when resuming.
General databases use persistent shared stores; project-specific protein sources
use `protein_source_store`. The copied final release is independent of backing
work files; intermediate linked outputs may not be.

## Supported scope and validation

The current prediction integration requires the standard nuclear genetic code.
Automatic taxonomy mappings are finite and require overrides when unmapped.
Pinned UTR ranking supports reverse/unstranded RNA; forward-stranded libraries
require `--utr_refine false` until an adapter is validated. Portcullis is
complementary QC, not a mandatory pre-BRAKER filter or AED component.

Regression coverage includes input contracts, native tool handoffs, optional
evidence modes, output consistency and cache behavior. Working code still needs
these tests when tools, settings or implementation change. See
[validation scope](docs/validation.md) for the distinction between synthetic
checks, software readiness and biological validation; it is not a guarantee for
every taxon or evidence combination.

[Original/current audit](docs/rework-audit.md) and [cleanup inventory](docs/cleanup.md)
record intentional changes, completed source consolidation and the retained validation layout.
