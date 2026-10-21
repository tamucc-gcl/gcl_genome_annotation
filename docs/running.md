# Running and resuming

## Execution environment

Use Linux with Nextflow DSL2, Java compatible with the installed Nextflow,
Conda and Singularity/Apptainer available through the site's Singularity command.
Nextflow 23.10.1 is the cluster-tested baseline; the repository does not establish
compatibility with every newer release. Python 3.11 is recommended for running
the test controllers; task interpreters and tools are configured separately.

`main.nf` is the normal entry. Invoke it from the consuming project, with an
explicit `-C /path/to/gcl_genome_annotation/nextflow.config`, `-profile slurm`,
`-w /persistent/project/work`, sheets and database arguments as shown in the
[README](../README.md). No manual stage approval/import is required.

The root config is currently Crest-oriented: `normal` queue, `/scratch`,
`/work/birdlab/.conda_builds` and `/work/birdlab/singularity_cache`, plus default
BUSCO/OMArk paths. Review/adapt these site settings before another deployment.
The `local` and `docker` profiles exist, but are not substitutes for validation
of the full large-genome workflow on another runtime. No personal sbatch script
or consuming-project inputs belong in the pipeline repository.

## Database paths

Pass store parents, not arbitrary internal data files, unless a parameter
explicitly requests a file. Mounts must be accessible to the driver and workers.

| Parameter | Existing layout / meaning |
|---|---|
| `taxdump_dir` | Fixed NCBI `names.dmp`, `nodes.dmp`, optionally `merged.dmp`; never downloaded by preflight |
| `orthodb_database` | General OrthoDB archive and validated-source store; see protein guide |
| `protein_source_store` | Project spike-in normalized-source store; defaults to launch-directory `protein_source_store/` |
| `dfam_database` | `dfam_db/dfam39_full.0.h5` through `.16.h5`, `RepeatPeps.lib`, derived subset/index stores |
| `diamond_database` | `swissprot.dmnd` for default functional search and default repeat-library filtering |
| `rfam_database` | `rfam_db/Rfam.cm`, `Rfam.clanin`, and `Rfam.cm.i1f/.i1i/.i1m/.i1p` |
| `busco_database` | `lineages/<resolved-lineage>/`, current default taxonomy mapping uses odb10 |
| `omark_database` | `<oma_db_name>.h5` (default `LUCA.h5`) and `taxa.sqlite` |
| `interproscan_database` | `interproscan-data-<interproscan_version>/data/` |
| `eggnog_database` | `eggnog_db/` with protein index, annotation and taxonomy databases |
| `kofamscan_database` | `kofamscan_db/ko_list` and `profiles/` |

Existing stores use native `storeDir` reuse, including fresh projects without
`-resume`. A configured release selector is not proof of the contents of an
already populated store. Keep stores immutable during runs; use distinct parents
for intentional database changes. Missing standard sources use the available
acquisition modules; inspect the relevant stage guide for validation requirements.
Main's default de novo repeat branch expects the general DIAMOND index to exist
unless `repeat_filter_database` supplies a different file. Functional DIAMOND
acquisition later in the workflow must not be relied upon to prepare that input.

## Evidence modes and switches

Assembly `rna_readsets` controls requested RNA; omit the RNA sheet only when no
assembly requests it. Use `use_protein_evidence=false` for RNA-only/genome-only
assemblies; otherwise taxonomy can select OrthoDB even without a spike-in.
`protein_fastas=none` alone does not disable OrthoDB. See [input precedence](foundation-preflight.md).

| Control | Effect |
|---|---|
| `repeat_library` | Supply a complete masking library instead of de novo discovery |
| `repeat_reclassify=false` | Skip optional Terrier classification |
| `utr_refine=false` | Keep original coding models; RNA evidence can still be assessed |
| `run_aed=false`, `run_junction_qc=false` | Disable the corresponding assessment; disable dependent selection rules too |
| `run_baseline_qc=false` | Omit unfiltered scientific QC; selected QC remains separate |
| `run_te_flagging=false` | Omit TE evidence work; cannot request TE removal |
| `run_functional_annotation=false` | Omit functional searches/decoration; release coding models remain available |
| `run_<tool>=false` | Disable BUSCO/OMArk, one functional tool or one independent track; supported row overrides are in the input guide |
| `run_finalization=false` | Omit final IDs, selected QC, release and unified report; also disable functional annotation |
| `run_model_selection=false` | Also requires finalization/function disabled |
| `run_model_assessment=false` | Also requires selection/finalization/function disabled |
| `braker_check_only=true` | Preparation/alignment and actual BRAKER software checks, without training/postprocessing/release |

`braker_check_only` is not a Nextflow stub run: upstream tools can still execute.
Omitted paths do not disable tools. Missing enabled results fail validation rather
than quietly dropping an assembly. Scientific tool/scheduler failures may stop
the run before a final summary; manifest-level isolation is not universal runtime isolation.

## Resources

Resource defaults are requests, not measured requirements or guaranteed memory
caps. They can be substantial: mapping reserves 16 CPUs/48 GB, BRAKER 48 CPUs/
200 GB, protein evidence alignment and native AED 64 CPUs/300 GB, UTR ranking
16 CPUs/64 GB, Portcullis 8 CPUs/48 GB, and centroAnno 64 CPUs/300 GB while using
8 native threads. See config for times and tool-specific settings.

The whole-node centroAnno reservation reflects Crest's memory-sharing history;
it is not a general requirement to run every tool on 64 threads. BRAKER tool
threads are capped at 48. HISAT2 defaults to `rna_no_temp_splicesite=true` after
runaway memory with within-run splice reuse; this changes sensitivity and is
recorded. Resources, task concurrency and disk needs must be assessed for the
actual dataset and scheduler. Initial environment/container/database setup can
dominate short acceptance runs; native AED against a large evidence set can
take more than a day.

## Resume and storage

Use normal `-resume` from the same launch directory, with the same work directory
and retained `.nextflow/` cache. A named earlier run can be selected explicitly
when needed; switching entry points/workflow aliases can change task hashes.
Preflight/planners may rerun deliberately to inspect current input state; this
does not mean all expensive tasks rerun. Compare process trace status and inputs.

| Change | Expected affected work |
|---|---|
| Report formatting | Reporting only |
| Final release export helper | Packaging, index and downstream report |
| Selection thresholds | Selection, final IDs/QC, affected searches/deliverables/reports; prediction and evidence searches reused |
| Track QC threshold | Track selection and dependent packaging/report; native track search reused |
| Supplied library strand | Strand-dependent QC/refinement/evidence; neutral mapping need not rerun |
| Reference, reads or mapping options | Affected indexing/mapping and scientific descendants |
| Final ID prefix | ID-dependent selected QC, functional queries/results and delivery |

Nextflow owns staging, scratch, transfers and cache. `storeDir` is a persistent
output/reuse destination, not the task working directory. Do not clean work or
source FASTQs while intermediate outputs still link to them. Release payloads
are copied and checksummed. Old published files are not automatically deleted;
use current manifests to export results. Never delete a Conda prefix while an
environment build or task is still using it.
