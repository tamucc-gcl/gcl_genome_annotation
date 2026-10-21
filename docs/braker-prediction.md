# BRAKER3 prediction

`main.nf` joins the accepted softmasked reference, per-assembly proteins and
separate RNA BAMs by explicit identity, checks the runtime and performs BRAKER3
training/prediction. Postprocessing, function, tracks, release and reporting are
connected downstream. Missing requested evidence blocks the job rather than
changing its mode silently.

| RNA | Proteins | Route | Evidence flags |
|---|---|---|---|
| Yes | Yes | ETP | `--bam` and `--prot_seq` |
| Yes | No | ET | `--bam` |
| No | Yes | EP | `--prot_seq` |
| No | No | ES | `--esmode` |

RNA associations and `use_protein_evidence` in the assembly sheet select modes.
BAMs use a stable readset order and remain separate; BRAKER handles transcript
assembly internally in ETP. Genome-only ES self-trains without a pretrained
species model. `braker_ab_initio=false` controls extra no-hint AUGUSTUS outputs,
not ES availability. The current integration requires nuclear genetic code 1.

## Defaults and validation

The configured image is `teambraker/braker3:v3.0.8`; defaults are 48 reserved
CPUs, at most 48 native threads, 200 GB, 240 hours and eight optimization rounds.
These requests are configurable, not a biological/runtime guarantee. Each job
gets a private writable AUGUSTUS config and deterministic assembly-specific
training name. Trained parameters and native work/logs are retained.

Input checks reconcile checksums, reference/BAM dictionaries, layout/groups,
readable indices and accepted alignment QC. Output checks require nonempty
predictions, reference bounds, valid feature IDs/parents and matching coding
transcript/protein/CDS IDs. Ready means integrity checks passed, not that every
gene is biologically correct. BRAKER UTR prediction and strand-separated coverage
BAM options are not enabled: resolved library strand is metadata used by later
RNA analysis, not passed as an unrelated BRAKER strand flag.

## Outputs and focused entry

`structural_annotation/README.md` and `prediction_status.json` describe modes and
readiness. Input/runtime records are under `braker_inputs/` and `braker_runtime/`.
`braker/<assembly>/` retains original GFF3, proteins, CDS, provenance, execution
logs, training parameters and `braker_run/`. Optional UTR work never overwrites
these predictions.

`validation/entrypoints/braker_prediction.nf` is a diagnostic entry that reads accepted reference/repeat,
protein and alignment reports under outdir, or `braker_inputs_root`. It stops at
prediction/reporting; it does not execute the rest of main. Unused evidence
reports can be absent. `braker_check_only=true` requests actual software checks
without training. Moving between entry points can change cache identity; normal
operation is the same `main.nf` invocation with `-resume`.

```bash
python gcl_genome_annotation/tests/validate_braker_cluster.py \
  --outdir testing/braker_contracts_01
```

The acceptance fixture runs actual upstream preparation/masking/mapping and
BRAKER software checks for all four evidence routes; it does not train models.
Real RNA+protein ETP training and downstream release completed. Full real-data
training of ET/EP/ES is not established by those software checks.
See [validation scope](validation.md) and
[BRAKER v3.0.8](https://github.com/Gaius-Augustus/BRAKER/blob/v3.0.8/README.md).
