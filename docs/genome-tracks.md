# Remaining genome annotation tracks

`main.nf` reconnects the existing tRNAscan-SE, parallel Infernal/Rfam, tidk and
centroAnno modules. These independent branches do not consume gene predictions,
functional-search results or RNA alignments. They run alongside those stages as
their references become available. BRAKER software-check-only runs omit them.

## Reference and settings

tRNAscan-SE, tidk and centroAnno receive the normalized, filtered **unmasked**
reference from `FILTER_ASSEMBLY`. Infernal retains the existing softmasked input
and whole-contig balanced chunks. All coordinates use the same retained contig
names and lengths as the coding annotation. Infernal uses the full assembly's
both-strand search space for E-values even when searches are chunked.

The existing per-assembly/global settings control `run_trnascan`, `run_infernal`,
`run_telomere` and `run_centromere`; their defaults remain enabled. tidk receives
each assembly's preflight-resolved `telomere_motif`, rather than a shared fallback.
Window/count/end-distance settings remain in the single `nextflow.config`.

The initial full-genome centroAnno pilot was SIGKILLed three times after reaching
203–221 GiB RSS despite a 16 GB request. Accounting did not report `OUT_OF_MEMORY`,
so this strongly suggests memory pressure without establishing the kill source.
The default now reserves a whole Crest normal node (`centromere_cpus=64`,
`centromere_memory='300 GB'`) while running at most `centromere_threads=8` native
threads. Tool threads are capped by the allocated CPUs, so small acceptance jobs
still use two. These are configurable resource settings, not a verified memory
cap or guarantee of completion. The eight-thread rerun also failed at the same
step. Inspection of the actual Conda artifact then identified a CLI mismatch:
the rework's `-G true` was based on an upstream tag that differs from this package.
The original `-x anno-asm` is restored; the corrected real-genome pilot passed.
Retain the conservative allocation until actual assembly-mode usage is measured.
The label terminates on failure rather than repeating identical memory-heavy
attempts. Other clusters can reduce the allocation through these same parameters.

`--rfam_database /work/birdlab/databases/rfam` is the **store parent**, containing
`rfam_db/Rfam.cm`, `Rfam.clanin` and `Rfam.cm.i1f`, `.i1i`, `.i1m`, `.i1p`.
The user confirmed this existing layout. Nextflow's native `storeDir` reuses it
across projects without requiring `-resume`. Missing stores use the original
download/cmpress module; incomplete indexes fail clearly before cmscan.
`rfam_version` remains the existing acquisition setting (`CURRENT` by default),
not a claim that an existing shared store is the current release. Changing that
setting does not replace a populated store; use a separate store for a new release.

## Reused modules and corrections

- tRNAscan-SE uses its supported eukaryotic `-E` mode; native isotype scanning
  remains enabled. The old unsupported `--eukaryotic`/`--isotype` switches were
  removed. Successful no-hit runs can emit empty optional tracks.
- tidk's TSV `window` column is the **end** coordinate. The converter checks
  complete ordered windows, clips the final window to the actual contig length,
  merges adjacent above-threshold windows and classifies terminal/interstitial
  regions against that length. It does not infer length from rounded windows.
- The existing qjh96 package is pinned to **centroAnno 1.0.2 build 0**. Inspection
  of its executable confirms `-x anno-asm` is its assembly mode, as in the original
  module; its default is `anno-sat-asm`. The rework's replacement with `-G true`
  was incorrect. The package recipe builds a developer-local source folder, not
  the upstream `v1.0.2` tag, so identical version numbers do not establish identical
  command-line behavior. Acceptance now checks native assembly-region scanning
  in the captured `<assembly>.centroanno.log`, not merely a zero exit and files.
  Native monomer/HOR CSV ends are zero-based **inclusive**. Both ends are converted
  to one-based inclusive GFF3, with declared row lengths checked. These are
  satellite/HOR candidates, not experimental evidence of active centromeres.
- The original Infernal converter is retained. Validated output preserves its
  Rfam family in `Name` and uses feature type `ncRNA`. tRNAscan and Rfam can report
  overlapping loci; no unvalidated deduplication or gene exclusion is introduced.

Native format references: [tidk window writer](https://github.com/tolkit/telomeric-identifier/blob/main/src/search.rs),
[tRNAscan-SE options](https://github.com/UCSC-LoweLab/tRNAscan-SE/blob/master/tRNAscan-SE.src),
[actual qjh96 package artifact](https://api.anaconda.org/download/qjh96/centroanno/1.0.2/linux-64/centroanno-1.0.2-0.conda),
[monomer writer](https://github.com/junhaiqi/centroAnno/blob/v1.0.2/include/moutput.h),
[HOR writer](https://github.com/junhaiqi/centroAnno/blob/v1.0.2/include/hor.h).

## Post-search QC defaults

All four tracks now pass through a separate configurable QC stage. The main
GFF3 filenames and workflow output channel contain the selected features;
validated unfiltered GFF3 files remain under `raw/`. Changing a filter setting
reuses native searches with `-resume`. The first run of this update needs one
tRNAscan-SE rerun to retain its score breakdown and secondary structures for
the native confidence filter. No other search command changes for this update.

| Track | Default selection | Configuration in `nextflow.config` |
|---|---|---|
| Telomere | At least 50% motif coverage in each native window; merged dense region within 10 kb of a contig end | `telomere_filter`, `telomere_filter_min_density=0.50`, `telomere_filter_max_end_distance=10000` |
| tRNA | Native EukHighConfidenceFilter `high confidence set`, with children following their parent | `trna_filter`; native score cutoffs listed below |
| Rfam ncRNA | Native inclusion `!`; remove lower-ranked native overlaps `=` | `ncrna_filter`, `ncrna_filter_remove_overlaps=true` |
| Centromeric repeats | Native identity at least 0.80; HORs additionally require native `high` confidence | `centromere_filter`, `centromere_filter_min_identity=0.80`, `centromere_filter_high_confidence_hor=true` |

Set any filter's enable parameter to `false` to use its validated raw features.
These switches do not disable the search itself. The existing `run_*` switches
still control whether a track runs at all.

Telomere density is motif count times motif length divided by actual window
length, including a clipped final window. Density is checked **before** adjacent
qualifying windows are merged. End distance is the smaller gap between the
merged region and either contig end, using the reference dictionary. Thus weak
flanks cannot dilute a dense array and the old coarse terminal/interstitial
label does not determine selection. The 50% and 10 kb defaults are transparent
starting heuristics, not established biological cutoffs. A fragmented assembly
can have interstitial repeats near artificial contig ends. Native tidk windows
and permissive raw regions remain available for inspection.

The tRNA stage runs the existing
[EukHighConfidenceFilter](https://github.com/UCSC-LoweLab/tRNAscan-SE/blob/master/EukHighConfidenceFilter.in),
then projects its classifications onto the native GFF3. It does not reimplement
the native copy-number or isotype filters. Default native cutoffs are
`trna_filter_cmscore=50`, `trna_filter_ssscore=10`, `trna_filter_isoscore=70`,
`trna_filter_tertiary_start=70` and `trna_filter_tertiary_max=95`.
This is a eukaryotic nuclear tRNA confidence workflow; the native filter excludes
its recognized mitochondrial sequence names (`chrM`, `M`, `chrMT`, `MT`).
Native classifications, structures and log are published under
`raw/trna_confidence/`. Low-confidence and pseudogene calls remain in raw output.

Rfam searches already use curated family-specific GA thresholds. Following the
[Rfam annotation guidance](https://docs.rfam.org/en/latest/genome-annotation.html),
post-search QC uses Infernal's inclusion and overlap annotations, without adding
an arbitrary global E-value threshold or custom overlap ranking. The existing
search annotates overlaps across models; this is not restricted to clan members.
tRNAscan and Rfam remain separate potentially overlapping tracks.

The pinned centroAnno build supplies identity and HOR confidence fields. QC
uses those existing fields without assuming human monomer length or one
centromere per contig. Identity 0.80 is a configurable heuristic; satellite/HOR
similarity does not demonstrate functional centromeres. No tool upgrade or
additional centromere caller is introduced.

## Published deliverables

Native outputs retain the existing `ncRNA/`, `telomeres/` and `centromere/` layouts.
Per-assembly validated deliverables are under `genome_tracks/<assembly_id>/`:

- `README.md`: enabled/disabled status, raw/selected counts, policies and links.
- `genome_tracks_status.json`: reference identity, settings and output checksums.
- `<assembly_id>.<track>.gff3`: selected features per enabled track, with namespaced IDs and
  consistent parent references.
- `raw/<assembly_id>.<track>.gff3`: validated unfiltered features.
- `<track>.status.json`: raw/selected checksums, counts and filter policy.
- `<track>_filter/decisions.tsv`, `rejected.gff3`, `policy.json`: selection audit.

Telomere decisions and rejected features are recorded per candidate window;
selected features are merged regions. Their counts therefore cannot be compared
by subtracting rejected windows from raw region counts. Other tracks report
feature-level decisions, including parent-dependent child decisions. Native
FASTA/BED outputs retain their raw meaning and are not selected sequence sets.

Validation blocks unknown contigs, out-of-bounds coordinates, invalid strands,
duplicate IDs, dangling parents and missing/duplicate enabled results. A completed
search with no features is distinct from a disabled track. Old published files can
remain after disabling a track; the current status manifest is authoritative.
Genome-track publication explicitly overwrites older reports/GFF3 files on resume,
including when cached enabled results replace a previously disabled report.
Native inputs stage under a separate path from normalized outputs, so normalization
cannot overwrite its source or report a checksum of modified input bytes.

These outputs are separate from the annotated coding GFF3. Combining accepted
tracks into the final release and unified reporting follows the native pilot;
coding identifiers, model structures and functional tables are unchanged here.
Reports are Markdown; the existing native tidk SVG remains available.

## Cluster acceptance

After syncing, run with the usual Nextflow/Conda modules loaded:

```bash
python gcl_genome_annotation/tests/validate_genome_tracks_cluster.py \
  --outdir testing/genome_track_filter_contracts_01
```

The runner uses actual tRNAscan-SE, tidk, centroAnno and a tiny locally built and
calibrated Infernal covariance model. It exercises two motifs/assemblies,
fewer contigs than requested chunks, enabled no-hit telomere/tRNA searches,
native tRNA confidence filtering, disabled tracks/filters, threshold changes,
native database-store reuse, unchanged resume and report-only changes.
It does not download Rfam or run repeats, mapping, BRAKER or functional searches.
Use a new test directory for this update because the fixtures changed. Then use
`--resume` with that same directory to recover a failed test.

Local Python and Nextflow parsing/routing checks are not a substitute for this
Linux-tool acceptance. The native search pilot passed previously; the new QC
stages still require this acceptance and a real-genome review. The existing
normal sbatch command needs no extra arguments: retain
`--rfam_database /work/birdlab/databases/rfam` and `-resume`.
