# Unfiltered annotation baseline

The first downstream increment reuses AGAT, BUSCO and OMArk from the original
pipeline. It assesses accepted BRAKER outputs without filtering or renaming them.
The integrated `main.nf` runs it after prediction validation when
`run_baseline_qc = true` (default). It is skipped in BRAKER software-check mode.
No manual checkpoint is required in a normal pipeline run.

AGAT assesses the complete GFF3, including alternative isoforms. BUSCO and OMArk
use one representative protein per coding gene: the longest predicted protein,
with ties broken by transcript ID. Terminal `*` is excluded from length ranking;
the sequence itself is retained unchanged, including partial/internal-stop models.
The selection uses explicit GFF3 Parent relationships, not transcript name suffixes,
and creates a mapping table. It does not remove genes or replace the full proteome.

BUSCO receives each assembly's resolved lineage; OMArk receives its effective
taxid. `run_busco` and `run_omark` are resolved per assembly. Disabled tools are
reported as `not_requested` and do not block joins or require dummy evidence files.
BUSCO retains its existing 5.3.2 tool version and odb10 lineage baseline; AGAT and
OMArk retain their existing environments. BUSCO's full results and OMArk's output
directory are preserved alongside human-readable summaries and version records.

The existing shared database modules use native `storeDir`; the downloads run
locally without scratch so the driver and transfer share filesystem metadata.
Existing `lineages/<lineage>/`, `<oma_db_name>.h5` and `taxa.sqlite` outputs skip
acquisition. The user confirmed `actinopterygii_odb10`, `LUCA.h5` and `taxa.sqlite`
are present on Crest. Scientific QC remains on Slurm. Database compatibility and
actual BUSCO/OMArk execution still require the real pilot.

All parameters, environments and process resources remain in `nextflow.config`.
No additional configuration files are needed.

## Cluster contracts

```bash
python gcl_genome_annotation/tests/validate_annotation_baseline_cluster.py \
  --outdir testing/annotation_baseline_contracts_01
```

This uses tiny synthetic validated predictions and actual AGAT. It checks
multi-assembly processing, representative selection, disabled optional metrics,
unchanged resume, addition of an assembly without rerunning existing AGAT tasks,
and blocking files whose checksums no longer match validation. It does not run
BRAKER, repeat discovery, alignments, BUSCO or OMArk. Interrupted acceptance can
be recovered with the same directory and `--resume`.

Local Python contracts passed; cluster acceptance remains pending.

## Real flounder pilot

After cluster contracts pass, change only the Nextflow script target in the
project's sbatch from `braker_prediction.nf` to `annotation_qc.nf`. Keep the
sample/RNA sheets, taxonomy/database arguments, existing output directory
`testing/genome_annotation` and `-resume`. This entry reads the completed
predictions and runs QC only; BRAKER, RNA preparation/alignment and repeats are
not invoked. `baseline_inputs_root` optionally points at a different accepted
output root; by default it is the current `outdir`.

Inspect `qc/baseline/README.md` and `qc/baseline/baseline_status.json`.
Per-assembly directories contain `selection/`, `agat/`, `busco/` and `omark/`.
The report includes links and BUSCO/OMArk text summaries. `ready` means requested
tools completed with nonempty outputs, not that a biological quality threshold
was met. Preserve published predictions and backing work directories when outputs
are links. The development entry checks ready prediction status and GFF3/protein
checksums and the genome's normalized contig identity before selecting representatives.

## AED and subsequent filtering

MAKER AED measures agreement with available transcript and protein evidence; a
protein-only score is an evidence-specific assessment, not a reproduction of
MAKER's combined-evidence score. The existing COMPUTE_AED module compares BRAKER
CDS coordinates against the best overlapping miniprot protein alignment. Retain
it explicitly as protein-based AED for historical comparison. The user clarified
that the final AED output must assess all supplied evidence, including RNA when
available. Its implementation and transcript assembly/polishing choices are deferred
until after baseline QC. Missing evidence is unavailable, never substituted with
1.0. Combined-evidence scoring requires separate validation before interpreting
the legacy 0.7 cutoff; protein-only scoring is not the final RNA-inclusive policy.

Sources: [MAKER-P paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3912085/),
[AGAT statistics implementation](https://github.com/NBISweden/AGAT/blob/master/bin/agat_sp_statistics.pl),
[OMArk usage](https://github.com/DessimozLab/OMArk#usage).
