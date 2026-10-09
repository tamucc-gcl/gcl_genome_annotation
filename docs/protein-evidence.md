# Prediction protein preparation (shared-store acceptance pending)

`protein_evidence.nf` runs the accepted input preflight, prepares requested protein
sources, combines them per assembly, writes reports and stops. It does not run RNA
processing, repeats, alignment, BRAKER, model filtering or functional annotation.
The legacy `main.nf` still runs its original biological workflow. This is the next implementation increment,
not a claim that the full rework is ready for production.

## Sources and validation

- OrthoDB: taxonomy-selected odb12 partition, with its published archive MD5 pinned
  in `assets/orthodb_sources.json`. The catalogue and URLs come from the
  [Greifswald OrthoDB 12 partitions](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/).
  An existing `<orthodb_database>/<clade>.fa.gz` can be staged and checked read-only.
  The uncompressed legacy `<clade>.fa` is not automatically reused because its
  filename alone does not establish its release/checksum. Otherwise download into
  the task work directory; no shared database is overwritten.
- Reference accession: pinned NCBI assembly accession from the resolved settings,
  acquired with `datasets download genome accession ... --include protein`.
  Require exactly `ncbi_dataset/data/<requested_accession>/protein.faa` in the
  package; do not select whichever protein file appears first. Preserve the ZIP,
  annotation metadata, command, download log and actual tool version. See
  [NCBI download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)
  and [package structure](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/).
- Local spike-ins: explicit resolved paths, including gzip FASTA; verify the staged
  file matches the planning checksum and preserve the original file.

Each source is acquired once per distinct identity, independently of the number
of consuming assemblies. Source-specific protein IDs avoid collisions between
OrthoDB and spike-ins. A TSV crosswalk preserves original IDs and full headers.
Normalization uppercases amino acids, removes one terminal `*`, and wraps records
at 60 characters. Reject empty records, duplicate IDs within a source, invalid
amino-acid characters and internal stops. Ambiguous amino acids are allowed.
This checks FASTA/sequence integrity; it does not establish biological suitability.

Combine complete sources in stable source-ID order. Keep all sequences, including
identical proteins from different sources: no biological deduplication or filtering
is introduced in this stage. Disabled proteins produce explicit `not_requested`
status, with no empty FASTA/sentinel.

## Failure and resume behavior

Expected acquisition and validation failures produce a failed source manifest and
error log, preserving raw/partial downloads while removing partial normalized
FASTA/map files. All assemblies requesting that source become `blocked`; unrelated
assemblies can finish this stage. The final assertion exits non-success after the
summary is written if requested protein preparation is incomplete. Unexpected
programming, scheduler, environment or infrastructure failures can terminate the
workflow before final reporting and are not converted to successful status data.

Permanent `storeDir` bundles reuse successful sources across fresh runs and projects, without `-resume`. Native `-resume` also reuses unchanged assembly combinations. Preflight
and source planning deliberately rerun to check current inputs. Stable value
descriptors, rather than newly generated descriptor-file timestamps, enter task
hashes. Only an assembly's requested source bundles enter its combination.
Reporting uses a separate helper so report edits do not invalidate source tasks.
Previously failed published sources receive a retry token so resume retries them
without rerunning successful sources. Preserve work and launch history, and use
the same output directory while retrying failures. Advanced checkpoint imports
are not implemented.

The environment pins Python 3.11.9 and conda-forge NCBI Datasets CLI 18.35.0,
available for linux-64 in the [package catalogue](https://anaconda.org/channels/conda-forge/packages/ncbi-datasets-cli/files).
Source tasks use 1 CPU, 32 GB and 48 hours; combination tasks use 1 CPU, 8 GB and
24 hours. Synthetic cluster execution passed; real-download resource needs remain
validation gates. Source acquisition and initial store promotion retain both the original archive and
normalized outputs; ensure adequate shared storage before acquisition. Subsequent runs link the stored sources rather than copying the raw database into each run or assembly task.

## Cluster acceptance before real downloads

From the project root, after syncing the pipeline:

```bash
module load nextflow/23.10.1
module load miniconda3
python gcl_genome_annotation/tests/validate_protein_cluster.py \
  --outdir testing/protein_contracts_01
```

Use a new empty shared directory. The runner executes all input and protein unit
tests, then tiny local-FASTA Nextflow fixtures for shared/combined sources, disabled
proteins, failed source isolation, unchanged resume, adding an assembly, and a fresh launch/work/output directory reusing the same store without resume. It
asserts stored sources are never reacquired, unchanged combination tasks are cached, and only the new assembly's
combination runs when it reuses an existing source. It downloads no biological
data; it does create the pinned worker environment. Real NCBI/OrthoDB network
acquisition remains a separate pilot check. Inspect `acceptance.json` and the
per-case console logs if any check fails. Local checks were limited to source
review, with execution performed on the cluster.

The user-supplied acceptance.json on 9 October 2026 reports all checks passed:
unit tests, shared sources, unchanged resume, adding an assembly using existing
sources, failed-source isolation, and disabled evidence. The blocked case correctly
exited 1 while its independent assembly completed. Real OrthoDB and NCBI downloads
and the real-data evidence bundle remain pending.

## Personal submission command

After fixture acceptance, change only these two parts of the existing Nextflow
command in the personal project sbatch:

```bash
-C "${PROJECT_ROOT}/gcl_genome_annotation/nextflow.config"
run "${PROJECT_ROOT}/gcl_genome_annotation/protein_evidence.nf"
```

Keep the existing sheet, RNA, taxid, taxonomy, database, output and resume arguments.
Keep the user's testing output and log directories. One explicit config entry is
used: the root `nextflow.config` holds every stage setting. No personal
paths or submission launcher are added to this repository.

The first real pilot requests the Vertebrata partition and GCF_024713975.1.
The provider lists Vertebrata as a 5.8 GB compressed download; that acquisition is
not a tiny smoke test. Some small OrthoDB partitions are documented by the provider
as requiring supplementary proteins; automatic selection does not prove that a
standalone partition is adequate for BRAKER.

## Outputs

Under the existing `--outdir`:

```text
preflight/README.md and resolved_inputs.json
protein_sources/<source_id>/
  request.json, source.json
  raw.fa, raw.fa.gz or raw.zip
  proteins.fa, id_map.tsv                 only when source is ready
  download.log/error.txt                  when applicable
protein_evidence/README.md
protein_evidence/protein_status.json
protein_evidence/assemblies/<assembly_id>/
  evidence.json
  proteins.fa, id_map.tsv                 only when evidence is ready
```

For the user's testing output root, the summary is
`testing/genome_annotation/protein_evidence/README.md`. A passing protein-stage
pilot means validated sources and evidence sets exist; it does not run annotation.

## Permanent database stores

General databases belong in their shared database directory. OrthoDB uses:

```text
<orthodb_database>/validated/odb12/normalization-v1/<source_id>/
    raw.fa.gz
    proteins.fa
    id_map.tsv
    request.json
    source.json
    STORE_READY.json
```

The identity includes the release, clade and pinned provider checksum, and the
store separates normalization versions. Two vertebrate projects requesting the
same partition reuse the same bundle. Planning checks completeness and raw/FASTA/map
checksums before allowing reuse. Failed acquisitions never enter the store;
incomplete or damaged stored bundles produce a clear error rather than silent reuse.
The successful-bundle process uses Nextflow `storeDir` and is skipped when its
required output already exists. Output `protein_sources/` entries are links to
those bundles. Assembly-specific combinations remain in each run's output.

NCBI spike-ins and local proteins use the project-level `protein_source_store`
parameter, default `${launchDir}/protein_source_store/normalization-v1/`.
Launch from the project root, or supply an absolute `--protein_source_store` when
launching elsewhere. The spike-in identity includes its exact accession; local
FASTA identities include their path and content checksum. These sources need not
be shared among unrelated projects.

The normal run downloads missing validated sources and stores them permanently.
One-time imports of earlier per-run downloads are not implemented. The old
unversioned `<clade>.fa` files are not relabeled as odb12 or overwritten.

The existing Dfam, Rfam, BUSCO, OMArk, InterProScan, eggNOG, DIAMOND and Kofam
acquisition modules already use their shared database parameters with `storeDir`.
All future general database acquisition in the rework must retain persistent,
version-specific shared storage; assembly-specific results stay in run outputs.

The earlier fixture acceptance remains recorded above. This shared-store extension
has been source-reviewed and needs a new cluster acceptance run, for example
`--outdir testing/protein_contracts_02`. The fixtures use a temporary project store
under that testing directory, with no real database downloads.
