# Reference and repeat preparation

`reference_repeats.nf` implements input preflight followed by reference normalization,
repeat discovery/curation and soft masking. It stops after the reference/repeat gate;
RNA alignment and gene prediction are not connected yet. Normal end-to-end execution
with native `-resume` remains the eventual default. No checkpoint imports are required.

The implementation reuses FILTER_ASSEMBLY, CURATE_DFAM, REPEAT_MODEL,
CURATE_REPEATS, CLASSIFY_REPEATS, RECLASSIFY_REPEATS, FILTER_REPEAT_LIBRARY,
COMBINE_REPEAT_LIBS and the existing parallel RepeatMasker workflow. Changes repair
integrity problems and adapt interfaces to per-assembly inputs. Nextflow owns staging,
scratch, publishing and cache reuse. Existing inputs are declared as path inputs.
The legacy main entry receives corresponding interface adaptations; it has not been
migrated to the complete new input contract or accepted as the reworked end-to-end run.

## Reference contract

- Resolve `min_contig_length` per assembly through preflight, default 1,000 bp;
  retain contigs at or above the threshold. Zero disables length exclusion.
- Validate DNA/IUPAC sequence, including gzip FASTA; reject missing/empty records,
  whitespace within sequences, invalid bases and gaps. Do not silently change lengths.
- Normalize unsafe characters in the first header token; preserve the full original
  header in a reversible TSV. Reject duplicate/colliding IDs, including excluded records.
- Uppercase the reference before masking. Write its FASTA index, ordered contig
  dictionary with sequence hashes, retained/excluded counts and bases, and excluded lengths.
- Accept masking only when contig IDs, order, lengths and uppercase sequence hashes
  match this normalized reference exactly. Soft masking may change letter case only.

Invalid sequence data blocks that assembly and appears in the stage report, while
independent valid assemblies can finish. Scheduler, environment or scientific tool
failures can still terminate the run; those are not converted into accepted results.
Incomplete reports produce a nonzero exit at the final assertion.

## Shared databases and discovery

For the flounder project, keep these existing arguments:

```bash
--dfam_database /work/birdlab/databases/dfam \
--diamond_database /work/birdlab/databases/diamond
```

The normal stage reads `dfam_db/dfam39_full.0.h5` through `.16.h5` directly and
uses `diamond/swissprot.dmnd`. It does not download replacement databases or copy
the full FamDB database into the project. Override the general protein index with
`--repeat_filter_database /path/to/custom.dmnd` when needed.

Taxon subsets are native storeDir outputs under
`dfam/subsets/dfam39/tetools-1.95/taxon_<taxonomy-hash>/`. The subset is checksummed;
FamDB release, partition inventory/sizes and native information output are recorded.
This does not claim complete checksums of the many large HDF5 partitions. Treat a
release directory as immutable. Initial extraction checks the declared release and
all 17 partitions. Reuse assumes the stored release has not been altered.

`dfam/RepeatPeps.lib` is reused through the original extraction process's storeDir.
Its derived index is stored under `dfam/indices/diamond-2.1.11/`, with a content-derived
name and source/index checksums. Assemblies and projects share that index.

Discovery retains native RepeatModeler consensus FASTA, Stockholm alignment, training
directories and logs, with an explicit seed (default 42). Record-safe curation removes
short, simple/low-complexity/satellite and duplicate consensuses, then retains optional
CD-HIT clustering. Classification and the two-protein-database filter remain:
RepeatPeps hits retain a family; sufficiently strong general-protein-only hits remove
it; no-hit families remain. Preserve excluded families, every filter decision, both
search outputs and source checksums. Search failures fail rather than becoming no hits.

Tool environments/resources and parameter defaults remain controlled by root
`nextflow.config`. [TEtools 1.95](https://github.com/Dfam-consortium/TETools/releases)
and [Dfam 3.9](https://www.dfam.org/releases/Dfam_3.9/families/FamDB/) are pinned.
The existing Terrier environment is now pinned to bio-terrier 0.4.0, torch 2.6.0 and
lightning 2.6.1: the old torch <2.6 constraint conflicts with the current package's
[requirements](https://github.com/rbturnbull/terrier/blob/main/pyproject.toml).
Terrier 0.4 uses updated training weights; this is recorded in the stage settings and
must be evaluated in the real pilot, rather than assumed equivalent to older models.
The [upstream release description](https://pypi.org/project/bio-terrier/) describes
that model change. `--repeat_reclassify false` skips this optional classification step.

## Masking and outputs

Batching balances whole contigs across at most `n_rm_chunks` batches. Every contig
occurs exactly once, and the actual batch count controls gathering. There are no
artificial tile boundaries. A large chromosome remains indivisible, so workloads
can be uneven; requested batch counts above the contig count do not create empty jobs.

Retain RepeatMasker's temporary short-ID mapping and restore original IDs safely,
including long names, internal double underscores and apparent coordinate suffixes.
Reconstruct FASTA in reference order and wrap only sequence lines. Linked repeat
fragments retain one linkage ID within each batch, then IDs are uniquely renumbered
across batches. GFF features have unique IDs, escaped attributes and linkage IDs.
Repeat interval coverage uses inclusive interval unions. Per-class unions can
overlap and must not be added together. The legacy `.tbl` category lengths remain
explicitly labelled alignment-fragment sums; the canonical unique coverage is in JSON.

Published locations beneath `outdir`:

| Location | Contents |
|---|---|
| `filtered_assemblies/<assembly>/` | Normalized FASTA, `.fai`, header map, excluded list, dictionary, statistics and status |
| `repeats/denovo/` | Discovery/training outputs, curation/classification/filter decisions and excluded families |
| `repeats/libraries/` | Combined masking libraries |
| `repeats/batches/` | Whole-contig batch accounting |
| `repeats/logs/` | Native masking logs |
| `repeats/maps/` | Temporary short-ID to reference-coordinate maps |
| `repeats/` | Soft-masked FASTA, repeat GFF/OUT, alignments, summaries and landscape |
| `reference_repeats/` | Stage README, per-assembly status/identities/union statistics; optional benchmark results |

Expensive outputs use native `-resume`; adding an assembly should reuse unchanged
assemblies. Changing only library filtering should preserve discovery outputs and
invalidate the filtering/masking descendants. Shared stores also work without
`-resume`. Keep the work directory and Nextflow cache for settings iteration.

## Cluster acceptance and real pilot

Implementation and static syntax checks are complete; runtime acceptance is pending.
The first cluster attempt passed Python contracts but stopped during module
compilation on an unescaped Python regex dollar anchor in the Groovy script.
That anchor now uses Python's `\Z` end-of-string form. A static interpolation audit
passed across all 13 imported process modules; a regression check now covers the
stage's include tree. No scientific tasks ran in that attempt, and corrected
Nextflow execution remains pending.
The second attempt reached the masking merge but failed because its TEtools
container provides `python3`, not `python`. Merge helper calls and the Dfam
extraction validator now use `python3`; other TEtools stage tasks were checked
for the same mismatch. Corrected cluster execution remains pending.
Run from the consuming project root after syncing the pipeline:

```bash
python gcl_genome_annotation/tests/validate_reference_repeats_cluster.py \
  --outdir testing/reference_repeat_contracts_01
```

Use an empty output directory. This runs all Python contract tests and actual
RepeatMasker/merge utilities with a tiny synthetic complete library. It checks
nonzero masking, exact comparison against a single unsplit batch, long-ID and
coordinate safety, per-assembly thresholds/zero, unchanged resume, adding an assembly,
invalid-reference isolation and invalid-library blocking. It requests four CPUs and
8 GB per masking job through root-config parameters. It writes `acceptance.json`,
per-case commands, console logs, Nextflow logs and published reports.

`--repeat_library /path/to/complete.lib` intentionally bypasses discovery, Dfam
extraction, classification and protein filtering. It is an explicit supplied-library
route, not fake discovery acceptance. `--repeat_mask_benchmark true` additionally
compares parallel masking with one whole-reference batch; leave it false for the
large real run to avoid duplicate masking.

After fixture acceptance, change only the pipeline entry in the personal sbatch
Nextflow command from `rna_preparation.nf` to `reference_repeats.nf`, retaining its
sample sheet, taxonomy/database arguments, testing outdir, profile and `-resume`.
Do not pass the synthetic `repeat_library` or benchmark option for the real pilot.
The flounder run exercises the complete discovery/classification/filter/masking
chain with existing shared databases and stops automatically before alignment.
Its runtime depends on discovery and genome masking, and can be substantially
longer than protein downloads or RNA preparation.

Real de novo tool execution, taxon extraction/store reuse, shared RepeatPeps index
reuse, downstream filter-change cache boundaries and repeat biology remain pending
until their cluster/pilot evidence is available. Passing this fixture alone does not
mark the full reference/repeat stage scientifically accepted.
