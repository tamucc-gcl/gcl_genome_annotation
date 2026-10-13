# Rework requirements: end-to-end use and baseline filtering

Accepted 9 October 2026. This repository-local addendum supersedes conflicting
checkpoint and default model-filter policies in the 25 September consolidated
build specification. The synced planning documents remain read-only references.

## Normal execution

The default is one uninterrupted end-to-end run. Users may rerun with Nextflow
`-resume` to change settings. Stage outputs and integrity validation are automatic;
no manual approval, file handoff or separate invocation is required between stages.
Optional checkpoint stop/import interfaces are deferred until the normal run and
native resume work reliably. Preserve useful stage manifests and expensive outputs
without implementing a second orchestration system.

Nextflow owns input staging, scratch directories, output transfers, publication,
and task caching. Pipeline helpers perform scientific validation and transformations;
they must not add file links/copies solely to pass existing inputs between tasks.
Local RNA validation and disabled transformations emit metadata while consumers
receive the existing FASTQs as declared path inputs.

Reuse the original pre-rework downstream modules and helpers wherever compatible.
Adapt their interfaces to the new input/evidence contracts and fix demonstrated
correctness issues; do not replace working repeat, prediction, filtering, or
functional annotation implementations simply because the orchestration is changing.

Filter-only changes must reuse prediction and unaffected evidence/search results.
Final selection, affected output projections, QC and reports may rerun. Changes to
scientific prediction inputs correctly invalidate prediction and its descendants.

## Default model filtering

Start evaluation with the existing biological filter policy and thresholds:

| Rule | Default threshold | Default action |
|---|---|---|
| Short single-exon CDS | <300 bp | Remove |
| Short gene span | <200 bp | Remove |
| All recorded introns short | All <10 bp | Remove |
| Available AED | >0.7 | Remove |
| TE evidence | Existing flag thresholds initially | Flag; retain |

Expose each rule independently, including threshold overrides, disable and
flag-only behavior. TE exclusion is opt-in. Do not add default exclusions for
partial CDS or internal-stop models; retain explicit flags and the planned
functional-search eligibility handling.

Correctness fixes still apply: missing AED is unavailable, never substituted with
1.0; do not sum CDS lengths across alternative isoforms; use explicit feature
relationships and consistent model/sequence sets. Document corrected transcript
and gene aggregation before implementing it. Equivalent thresholds do not imply
bit-for-bit reproduction of old buggy decisions. RNA-inclusive AED is a changed
score and requires validation before interpreting 0.7 as directly comparable.

Preserve candidates, excluded models and all filter reasons. Report retained and
removed counts, length/exon distributions, evidence availability/support and
enabled completeness metrics so the user can assess the baseline and alternatives.

## Implementation status

This records accepted behavior, not completed filter/checkpoint implementation.
The preflight currently resolves settings and file presence only; the legacy
annotation entry point has not been migrated. Real pilot inputs and personal
submission scripts are maintained separately in the consuming project repository.

Step 1 input-contract acceptance passed on 9 October 2026 under Nextflow 23.10.1:
24 Python tests plus ready, blocked and malformed Nextflow fixtures passed.
This supplements the successful real-data preflight. Read contents, source
acquisition and biological stages remain outside this acceptance scope.

Protein-source preparation is now implemented as a separate entry point, with
synthetic cluster acceptance passed on 9 October 2026 (user-supplied acceptance.json). See [protein evidence](protein-evidence.md) for tests
and the real-download pilot. The real protein pilot and stored-source rerun passed on 9 October 2026. Next: RNA
preparation, reference/repeat integrity, four-route prediction, model/evidence assessment, final filtering
and functional/release integration. The native resume dependency checks remain
required before declaring settings iteration reliable.

All pipeline parameter defaults, process environments/resources, execution profiles and reports now live in the root nextflow.config. Stage-specific config files have been removed. Species-dependent defaults are unset for taxonomy resolution; the preflight contig default is 1,000 bp. The legacy main.nf entry now reads its defaults from that same config; legacy runs must explicitly provide species settings because that entry does not yet consume the taxonomy-resolved manifest. The shared-config protein fixture suite passed on the cluster, including native resume and adding an assembly. Real OrthoDB and exact-accession spike-in acquisition and stored-source rerun passed; broader fresh-project reuse fixtures remain unconfirmed.

General databases must persist in shared, version-specific storeDir locations and be reused by fresh runs and different projects when the required validated source exists. Protein preparation now stores OrthoDB bundles under orthodb_database/validated/odb12/normalization-v1; exact-accession spike-ins/local proteins use a project-level protein_source_store. Failed sources are never promoted. Per user direction, do not implement one-time import of earlier per-run downloads. This extension requires a fresh cluster acceptance run, including reuse without -resume.

Protein real-data acceptance: user confirmed the first run succeeded and supplied the successful rerun logs annotation-1515087.out and .nextflow_annotation_1515087.log. Both protein sources loaded from stores; COMBINE_PROTEIN_EVIDENCE, REPORT_PROTEIN_EVIDENCE and ASSERT_PROTEIN_EVIDENCE loaded from native cache. No sources were reacquired. This confirms the real pilot and same-project stored reuse; it does not claim the failed store fixture suite was rerun successfully.

RNA acquisition/preparation is implemented with local and accession paired/single routes, full FASTQ integrity checks, raw/processed FastQC, conservative fastp trimming, pair-preserving deterministic subsets and MultiQC. All settings remain in the root config. Automatic strandedness is pending reference alignment, never silently unstranded. RNA Python syntax checks passed; unit/Nextflow cluster acceptance and the real RNA pilot remain pending. See docs/rna-preparation.md.

RNA first cluster attempt: unit_tests passed; acquisition environment solve failed before actual RNA process acceptance. Corrected Conda channel order/overrides and conditional acquisition environment (Python-only local reads; SRA tools only for accession reads). Tool pins are retained; corrected cluster acceptance remains pending.

RNA synthetic cluster acceptance passed on 9 October 2026 in testing/rna_contracts_02: unit tests, fresh paired/single/shared RNA processing, unchanged resume, subset-only rerun frontier, shared-readset failure isolation with expected exit 1, and no-RNA routing. Conda solve/runtime corrections are validated. Next is the real paired-end RNA pilot; accession network acquisition and reference-specific strand inference remain unvalidated.

Subsequent native-staging cleanup: local/passthrough RNA steps emit metadata only;
real reads are declared path inputs. FastQC naming uses native staging; reporting
stages manifests rather than read bundles. Protein store/publication handoffs use
native input/output staging instead of helper copytree/symlink operations. Actual
processing inherits configured scratch; shared-path planners remain outside scratch.
Python syntax checks passed. Updated RNA and protein handoffs await cluster validation;
prior acceptance applies to the preceding implementation. Original downstream
modules remain unchanged and are the starting point for subsequent integration.

The rna_contracts_03 cluster attempt stopped at module parsing: dynamic stageAs
option strings referenced rid before task inputs were bound. Replaced those
options with native dynamic path declarations for RNA QC and protein handoffs.
No RNA tasks ran in that attempt; corrected Nextflow execution remains pending.

RNA native-staging cleanup acceptance passed on 9 October 2026 in
testing/rna_contracts_04 (user-supplied acceptance.json): unit tests, fresh
paired/single/shared processing, unchanged resume, subset-only reruns, metadata-only
passthrough fresh/resume, expected blocked exit 1, and no-RNA routing. RNA staging
and scratch changes are accepted. Next is the real RNA pilot. The separately changed
protein store/publication handoffs, real accession acquisition, and reference strand
inference remain pending validation.

Real RNA local pilot: user supplied ready manifests for both libraries and confirmed
the combined MultiQC report displays correctly. It contains eight raw/processed
FastQC samples and fastp results. Split reporting now gives separate raw and filtered
MultiQC pages; filtered includes fastp when trimming runs. Report-content acceptance
checks now verify expected sample identities and stage separation. These reporting
changes were confirmed working in the real local RNA pilot by the user;
acquisition/trimming/subsetting code is unchanged. Updated synthetic report-content
checks have not been rerun on the cluster.

Reference/repeat preparation is implemented in reference_repeats.nf, using existing
discovery/classification/filter/masking modules with corrected per-assembly reference
contracts, record-safe repeat filtering, whole-contig batches, exact reconstruction,
linkage/GFF integrity and overlap-aware coverage. Existing Dfam 3.9 partitions,
SwissProt DIAMOND and RepeatPeps are reused; taxon subsets and derived RepeatPeps
indices use shared stores. Tool pins, including the repaired Terrier environment,
are described in docs/reference-repeats.md. Cluster fixture acceptance passed on
9 October 2026 under Nextflow 23.10.1 in testing/reference_repeat_contracts_03
(user-supplied acceptance.json): Python contracts, actual custom-library masking
with unsplit comparison, unchanged resume, adding an assembly, independent valid
reference completion despite an invalid reference, and invalid-library blocking.
The two blocked cases returned the expected exit 1. This accepts reference integrity
and the custom-library masking route. Full real-data discovery, classification,
protein filtering, shared repeat-database stores and repeat biology remain pending
the flounder pilot; fixture acceptance does not imply those branches were exercised.

The first real repeat pilot stopped at CURATE_DFAM because species taxid 141886
was not present in the Dfam taxonomy. The automatic taxa_name default now selects
the nearest ancestry group from the explicit catalogue (Actinopterygii for this
flounder), retaining the effective species taxid. Explicit global/row choices still
override it. RepeatModeler was gated by the failed extraction. Resume preserves
completed work when task inputs are unchanged; corrected extraction and real
discovery acceptance remain pending on the cluster.

12 October 2026: the user supplied completed real flounder reference/repeat
reports. Reference identity and all 679,059,826 bp across 31 contigs were
preserved. Unique masked coverage was 47.55%; 1,364 de novo families survived
protein filtering, with 61 gene-derived candidates removed and 307 TE-protected
candidates retained. The reports support continuing with this baseline;
near-completely masked small contigs still need downstream RNA/protein support
review. Fixture acceptance alone did not establish this biological result.

Next implementation is RNA alignment, followed by migration of accepted stages
into main.nf and extension of the existing BRAKER3 prediction module. The user
confirmed keeping BRAKER3; a conversational reference to BRAKER4 was incorrect.
BRAKER3 consumes multiple mapped BAMs and runs transcript assembly internally.
The alignment stage therefore produces separate sorted BAMs per readset/reference,
with no pipeline BAM merge or standalone StringTie run. Library associations
remain metadata for QC and strand consistency. The development entry reuses the
accepted RNA/reference outputs and does not run repeats again. This does not
introduce a checkpoint/import system into the eventual end-to-end workflow.
Python checks pass; actual Nextflow/HISAT2 cluster acceptance remains pending.

12 October 2026: RNA alignment cluster acceptance passed in
testing/rna_alignment_contracts_01 (user-supplied acceptance.json). After the
interrupted Python/pysam environment was rebuilt, recovery, unchanged resume,
strand-only changes, an added reference, expected blocked exit 1 and no-RNA
routing all passed using actual HISAT2 alignment. The environment versions were
compatible; the earlier build stall's cause remains unconfirmed. This accepts
the synthetic alignment/BAM handoff contracts, not real flounder mapping quality.
Next is the real alignment pilot using accepted prepared reads and references.

Real alignment retries: the first attempt lost both jobs to NODE_FAIL on
crest-c003. The next attempt placed both jobs on crest-c019 with 48 GB requests;
Slurm reported approximately 156 and 212 GiB peaks, one FAILED job and one
CANCELLED job, without an explicit OUT_OF_MEMORY state. The cause of the SAM
stream truncation remains unconfirmed. HISAT2 stderr now also reaches the native
task error log. The user reports unreliable exclusive-node behaviour on Crest
and confirmed that normal nodes have 64 CPUs and 350952 MB configured memory.
Real alignment now reserves 64 CPUs and requests 300 GB per task, while keeping
14 mapper threads. This uses the full CPU allocation to avoid co-locating the
alignment jobs; no exclusive-node flags are used. Synthetic acceptance uses
four CPUs, two mapper threads and 8 GB. The revised real pilot remains pending.

12 October 2026: a further full-node attempt reported HISAT2 killed with signal
9 and roughly 375 GiB task peak. A mapper-only two-thread IBC1 diagnostic still
grew to approximately 168 GiB before cancellation. With only
`--no-temp-splicesite` added, the full 58,267,066-pair diagnostic completed in
29m35s with 1,123,308 KiB peak RSS, 76.63% overall alignment and exit 0. Disable
within-run splice-site reuse by default through `rna_no_temp_splicesite = true`;
keep spliced mapping, DTA and XS tags. The setting remains configurable because
reuse can help short-anchor reads. The comparison strongly implicates reuse in
memory growth; it does not prove a particular software defect or equal biological
sensitivity. At the user's request, return alignment resources to 16 reserved
CPUs, 14 mapper threads and 48 GB memory. Full BAM/QC validation at these
resources remains pending; the diagnostic used two threads and discarded SAM.

RNA alignment reporting now generates MultiQC 1.25.1 from native staged HISAT2,
samtools flagstat and samtools stats text logs. Assembly/readset labels prevent
shared log filenames or reuse of a readset across assemblies from collapsing
samples. The generated README links to multiqc/multiqc_report.html; no-RNA
inputs receive a clear no-alignment page. Reporting has a separate helper and
environment so this change does not invalidate completed scientific tasks.
Local validation with actual MultiQC retained both assembly identities in all
three datasets; 70 Python contracts passed. Cluster integration checks now also
require the linked report and complete labelled MultiQC datasets.

12 October 2026: the real alignment outputs and split preparation/alignment
MultiQC reports were reviewed and accepted. IBC1/IBC2 primary mapping was
76.63%/80.92%, with tens of millions of spliced primary reads. Strand inference
was confidently reverse for both libraries. BAM dictionaries, indices, mate
counts and reference identities passed the existing reports. This accepts the
evidence for a BRAKER3 pilot, not gene-model quality.

The next implementation connects protein/RNA/reference preparation, repeats and
alignment in `main.nf`, then extends the existing BRAKER module for RNA+protein
prediction using multiple separate BAMs. A development `braker_prediction.nf`
reuses completed stage outputs to avoid rerunning the long real repeat/alignment
pilot. Runtime checks run automatically before prediction; optional check-only
execution is for diagnosis. Root config remains the sole configuration file.
BRAKER3 v3.0.8 is pinned, private AUGUSTUS training is retained, and the existing
eight optimization rounds remain the baseline. Library protocols are metadata;
BRAKER's strand-separated UTR BAM flag is not used. The current connected milestone
stops after prediction integrity reporting. Scientific downstream modules remain
available for the next filtering/annotation steps. Python contract validation
passes locally; connected cluster acceptance and real BRAKER training are pending.

BRAKER handoff cluster acceptance subsequently passed on Crest, Nextflow 23.10.1,
in `testing/braker_contracts_01` (user-supplied acceptance.json). Python contracts,
integrated recovery and unchanged resume, and prediction-only recovery and
unchanged resume all passed. The fixture used actual preparation/alignment tools
and BRAKER3 software checks; it did not train or predict genes. Corrections
separated the protein store output from its staged input, preserved native reference
paths through the connected planner, and distinguished first standalone recovery
from unchanged resume in the validator. Proceed to the real-data prediction-only
pilot using the accepted published stage outputs. Real training and biological
prediction quality remain unvalidated.

The user reaffirmed the original requirement to support every RNA/protein evidence
combination. The initial RNA+protein-only restriction was incomplete and is removed.
BRAKER now selects ETP, ET, EP or explicit `--esmode` per assembly. Optional evidence
files are staged only when requested; missing requested evidence blocks rather
than downgrades a route. The integrated join now retains RNA-free assemblies.
Prediction-only genome runs require only accepted reference/repeat outputs, with
no placeholder RNA/protein reports. Genome-only training uses BRAKER's documented
GeneMark-ES/AUGUSTUS route. The additional AUGUSTUS output flag remains configurable
and defaults false; ES prediction is intrinsically ab initio. Expanded Python tests
pass locally. The cluster validator now covers mixed evidence modes and a completely
genome-only workflow with unchanged resume. Its updated fixtures require a new
output directory; expanded cluster and real training acceptance remain pending.
