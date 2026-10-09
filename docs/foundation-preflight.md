# Foundation implementation: inputs and taxonomy

This is the first implementation increment of the agreed annotation refactor. `preflight.nf` resolves inputs and settings only. It does not acquire reads or launch repeat masking, alignment, BRAKER, functional annotation or webshare. The legacy `main.nf` has not yet been migrated to consume this manifest. Do not pass the new wide sheets to the legacy entry point expecting the new behavior.

## Precedence

1. For taxid, use a nonblank assembly-row value, otherwise the global `params.taxid`.
2. For each species-related setting, use an explicit row value, then an explicit global value, then derive it from the effective taxid.
3. An explicit `auto` at row level requests derivation even when a global override exists. Blank inherits. `false`, `0` and `none` remain explicit values.

Thus one global taxid is enough for a single-species batch, with taxid omitted from the assembly sheet. Different assemblies can override it. Per-setting overrides such as a global BUSCO lineage intentionally remain authoritative even for a row with a different taxid; use row `auto` to request derivation for that row.

`ncbi_speciesID` is accepted as a legacy alias for taxid, and `assembly_path` as an alias for assembly_fasta. Conflicting aliases are errors.

Derived settings:

| Setting | Source |
|---|---|
| Scientific name | NCBI scientific name for the effective node |
| OrthoDB partition | Closest matching ancestor in `assets/taxonomy_rules.json`; no automatic Eukaryota fallback |
| BUSCO lineage | Closest mapped ancestor; initial mapping targets odb10, matching the legacy pipeline generation |
| Dfam selection (`taxa_name`) | Effective taxonomy ID |
| Nuclear genetic code | Nuclear field in nodes.dmp, following its inheritance flag; not the mitochondrial code |
| Telomere motif | Conservative curated mapping; initial automatic coverage is vertebrates only |
| OMArk organism identity | Effective taxid |

For an unmapped telomere motif, provide an override or disable telomere analysis. Do not silently apply the fish motif to every organism. The mapping can be expanded with supported clade-specific rules and tests. BUSCO/OrthoDB dataset availability and tool compatibility remain cluster validation gates. Taxonomy does not identify the correct close-relative spike-in proteome: that remains explicit.

eggNOG's taxonomic restriction is intentionally unrestricted (`none`) unless overridden. A restriction is not equivalent to the organism's identity: automatically restricting a non-model genome to its own species could defeat the intended homolog search. All resolved settings and origins are saved for review.

## Cluster invocation

Run from your project submission directory with its paths available to compute nodes. Supply an existing fixed NCBI taxdump directory containing names.dmp and nodes.dmp, optionally merged.dmp. The resolver records checksums and never downloads or updates the taxonomy snapshot.

```bash
nextflow -C gcl_genome_annotation/conf/preflight.config \
  run gcl_genome_annotation/preflight.nf \
  -profile slurm \
  --sample_sheet data/annotation_assemblies.csv \
  --rna_sheet data/annotation_rna.csv \
  --taxid 373251 \
  --taxdump_dir /path/to/pinned/taxdump \
  --outdir annotation_preflight
```

`-C` is deliberate: it keeps legacy fish defaults and unrelated process configuration out of this new entry point. Add cluster-specific queue/account/conda cache settings through a separate configuration as needed. The preflight configuration uses SLURM, one CPU, 8 GB and a pinned Python interpreter; it needs no bioinformatics packages. Input FASTA/FASTQ paths are checked on the shared cluster filesystem, not copied into the preflight task. Relative paths inside a sheet resolve against that sheet, while global protein paths resolve against the launch directory. Inspect the resolved absolute paths.

Use global settings in a Nextflow params file or configuration exactly as on the CLI; the entry forwards the resolved parameter map. A row override always wins for that field. RNA rows specify local FASTQ or one SRR/ERR/DRR run accession; accession downloads and actual FASTQ-content/layout validation come in the RNA increment.

Outputs under `annotation_preflight/preflight/`:

- `README.md`: assembly readiness/mode/reason summary, plus failed RNA inputs.
- `resolved_inputs.json`: resolved settings, origins, RNA associations, taxonomy lineage, snapshot/rules/sheet checksums and pending-validation labels.

A separate assertion task returns non-success after the report task when any requested input is blocked. Ready assemblies remain listed, but this entry launches no annotation branches. This is not yet the complete runtime failure-isolation implementation. Global malformed/ambiguous tables produce a failed preflight rather than a silently reduced batch.

## Cluster checks before using the contract downstream

```bash
cd gcl_genome_annotation
python -m unittest discover -s tests -p 'test_annotation_inputs.py'
```

These tests use tiny synthetic fixtures for precedence, quoted CSV, all four routes, shared RNA, blocked assemblies, source conflicts and nuclear genetic-code selection. Then run the Nextflow preflight on a small real cluster sheet and check global taxid inheritance, row taxid override, row `auto` against global overrides, and missing-taxonomy failures.

The user has successfully run the real-data preflight on Nextflow 23.10.1 with
an assembly, RNA inputs, taxonomy and spike-in accession, reporting `ready` and
`rna_protein`. A malformed RNA CSV also produced the expected published failure
report and non-success assertion. This does not validate biological stages.

Complete the remaining foundation acceptance checks from the project root:

```bash
module load nextflow/23.10.1
module load miniconda3
python gcl_genome_annotation/tests/validate_preflight_cluster.py \
  --outdir testing/preflight_contracts_01
```

Use a new empty output directory on the shared cluster filesystem. The runner
first executes the Python contract tests, then submits three tiny Nextflow
preflights: a ready batch containing all four routes and shared RNA; failed shared
RNA blocking only its consumers in the manifest; and malformed CSV with a global
failure. Non-success Nextflow exits in the latter two cases are expected and
checked. Reports, logs, synthetic inputs and `acceptance.json` are retained.
It launches no reads processing or annotation tools and leaves the real pilot
inputs, output directory, submission script and launch history untouched.

Run this controller command inside a cluster allocation if site policy prohibits
controllers on login nodes. Worker tasks use the existing SLURM preflight config.
The suite is stdlib-only; the controller needs Python 3.8 or newer. The newly
expanded tests and runner have received source review only and have not been
executed locally; their cluster result is the remaining step-1 acceptance gate.

## Next increment

See the [9 October rework addendum](rework-decisions.md) for the accepted default
filtering policy and simplified checkpoint scope. Real input sheets, cluster paths
and personal submission scripts belong in the consuming project repository.

After acceptance passes, implement protein-source preparation, then RNA
preparation, then reference/repeat integrity, following the revised implementation
order. Keep species settings out of unrelated read preprocessing and pass only
stage-relevant values to expensive tasks. The preflight entry remains separate
until the annotation stages are migrated.
