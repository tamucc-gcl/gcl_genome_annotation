# BRAKER3 prediction

This stage joins the accepted softmasked reference, combined proteins and separate
coordinate-sorted RNA readset BAMs by assembly when requested. It selects the
BRAKER training route from the evidence available for each assembly. There is no pipeline BAM
merge or standalone StringTie stage. The current milestone stops after validated
coding-gene predictions; downstream model filters, functional annotation and
biological evaluation are still to be connected.

## First cluster acceptance

RNA+protein acceptance passed on Crest with Nextflow 23.10.1 in
`testing/braker_contracts_01`: Python contracts, connected recovery/resume and
prediction-only recovery/resume all passed. This confirms the handoff and actual
BRAKER3 software checks, not training, predictions or biological model quality.
Expanded four-mode and genome-only acceptance passed on Crest on 12 October
2026 in `testing/braker_contracts_04`. All eight fresh/resume runs passed across
the connected and prediction-only entry points, including genome-only inputs
without optional evidence reports. These runs checked actual BRAKER3 software
readiness, not training or prediction.

From the consuming project directory, using Nextflow 23.10.1 and the same
Slurm/Conda/Singularity setup as earlier stage validators:

```bash
python gcl_genome_annotation/tests/validate_braker_cluster.py \
  --outdir testing/braker_contracts_04
```

The fixture runs the connected `main.nf` with tiny local proteins, synthetic reads
and a supplied repeat library. It runs actual preparation/QC, RepeatMasker and
HISAT2, validates the separate BAM handoff, and calls the real BRAKER3
`--checkSoftware` command. It then checks unchanged resume and the prediction-only
entry. The mixed sample sheet includes all four modes and checks each mode's
evidence flags. Additional genome-only runs have no RNA sheet; the standalone
genome-only input root has no RNA/protein status reports or files. Unchanged resume
must cache the handoff and runtime checks in each case. Use a new directory for this
expanded fixture; the old two-assembly fixtures cannot be resumed with it.
This test downloads the BRAKER container if needed, but does not run de novo
repeat discovery or gene training. Expect container/environment setup and queue
time to dominate. `acceptance.json` records success/failure and console/log files
remain in the fixture directory. An interrupted test can be recovered with the
same directory and `--resume`.

Python tests cover evidence mismatch blocking, unsupported modes, command
construction, separate BAMs, output coordinate/parent checks and matching predicted
protein/CDS identifiers. They do not demonstrate biological prediction quality.

## Real-data pilot

The user-supplied report on 13 October 2026 confirms that the real
`Sfl-n_primary` pilot completed in RNA+protein (ETP) mode and passed prediction
integrity validation. It produced 21,528 genes and 45,399 coding transcripts,
with 45,399 matching protein and CDS records. Both separate paired-end readsets
(`IBC1_S123_L006`, `IBC2_S124_L006`) were included with reverse strand metadata;
the input protein evidence contained 19,445,987 sequences. Validation covered
nonempty outputs, reference bounds, feature IDs/parents and matching sequence IDs.
These are unfiltered BRAKER outputs; biological quality and downstream filtering
remain to be assessed. Actual training of the other three evidence modes has not
been validated by this pilot.

Keep the existing sbatch setup and Nextflow arguments. Change the script argument
from `rna_alignment.nf` to `braker_prediction.nf`, retaining `-resume` and the same
`--outdir testing/genome_annotation`. No new config is required.

The prediction-only entry resolves the current sample sheets and loads:

- `reference_repeats/reference_repeat_status.json` and `repeats/<assembly>.masked.fasta`;
- `filtered_assemblies/<assembly>/<assembly>.contigs.json`;
- `protein_evidence/protein_status.json` and `protein_evidence/assemblies/<assembly>/proteins.fa` when proteins are requested;
- `rna_alignment/rna_alignment_status.json` and its listed BAMs/CSI indices when RNA is requested.

Optional status reports may be absent when their evidence type is unused by every
requested assembly. Requested evidence that is missing or invalid is blocked;
the pipeline does not silently downgrade the selected mode.

It does not rerun protein acquisition, repeats, preparation or alignment. It
checks reference/protein checksums, normalized reference identity, BAM dictionaries,
coordinate sorting, read groups and readable indices, retaining the previously
accepted streaming alignment QC. A BRAKER software check precedes training
automatically. Normal execution proceeds through full BRAKER training and
prediction, output validation and reporting, then exits; there is no need to stop it
at a checkpoint. `--braker_check_only true` is an optional diagnostic that exits
after the handoff/runtime checks.

If inputs are in a different output root, use `--braker_inputs_root <previous-root>`.
Retain the published inputs and their backing Nextflow work directories when they
are links. Native Nextflow manages staging, scratch and stage-out.

## Defaults and evidence handling

All settings are in root `nextflow.config`: pinned
`docker://teambraker/braker3:v3.0.8`, 48 reserved CPUs, up to 48 tool threads,
200 GB requested memory, 240 hours and the existing eight AUGUSTUS optimization
rounds. Memory/runtime need assessment on the real pilot; these are requests,
not measured requirements. Each task receives a private writable AUGUSTUS config
and assembly-specific training name. Trained parameters and BRAKER working outputs
are retained with the predictions.

BAMs are passed as a comma-separated list in stable readset order. Strand metadata
is retained, and the accepted HISAT2 splice-derived XS tags are used. BRAKER's
`--stranded` option refers to strand-separated coverage BAMs for UTR prediction;
it does not encode a reverse-stranded library protocol. UTR prediction is disabled.
The standard genetic code remains required. Per-assembly routes are:

| RNA readsets | Protein evidence enabled | Route | BRAKER evidence arguments |
|---|---|---|---|
| Present | Yes | ETP | `--bam` and `--prot_seq` |
| Present | No | ET | `--bam` |
| Absent | Yes | EP | `--prot_seq` |
| Absent | No | ES | `--esmode` |

Use `rna_readsets` and `use_protein_evidence` in the assembly sample sheet (or the
existing global protein setting) to select evidence. Enabled proteins must have a
valid configured source. No new mode config files or manual route switches are
needed. Genome-only ES trains GeneMark-ES on the assembly and trains AUGUSTUS
from those predictions; it needs no pretrained species setting.
`braker_ab_initio` remains the optional additional AUGUSTUS output flag for
evidence-supported runs, default false. Genome-only predictions are inherently
ab initio regardless of that additional-output flag. The
[BRAKER3 v3.0.8 source](https://github.com/Gaius-Augustus/BRAKER/blob/v3.0.8/scripts/braker.pl)
documents these options and the incomplete nonstandard-code support. The
[official genome-only example](https://github.com/Gaius-Augustus/BRAKER/blob/v3.0.8/README.md#testing-braker-with-genome-sequence)
uses `--esmode`.

## Results and subsequent runs

`structural_annotation/README.md` links to prediction files, runtime/input reports,
logs and RNA alignment MultiQC. `prediction_status.json` distinguishes software
checks from completed predictions. Per-assembly results are in
`structural_annotation/braker/<assembly>/`: GFF3, proteins, CDS, execution provenance,
trained AUGUSTUS parameters and `braker_run/` working outputs.

Integrity checks require nonempty predictions, valid reference bounds and feature
parents, and matching coding-transcript/protein/CDS IDs. Gene counts, completeness,
RNA/protein support, TE contamination and the later filter baseline still need
biological assessment.

For subsequent end-to-end runs, `main.nf` connects the same processes with native
output channels. RNA alignment can run alongside repeat preparation once the
normalized reference and prepared reads are ready; BRAKER waits for both evidence
branches and masking. Moving between standalone entry points and `main.nf` does
not guarantee reuse of their separate Nextflow cache sessions. Use the
prediction-only entry for this first pilot to protect the accepted repeat/alignment
work; normal `-resume` within either entry reuses unchanged tasks.
