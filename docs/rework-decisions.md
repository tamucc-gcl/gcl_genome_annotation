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

12 October 2026: expanded BRAKER acceptance passed on Crest with Nextflow 23.10.1
in `testing/braker_contracts_04` (user-supplied acceptance.json). Python tests and
all eight fresh/resume runs passed: mixed ETP/ET/EP/ES assemblies through both
entry points, plus genome-only connected and prediction-only runs without optional
evidence reports. The fresh protein-store handoff passed with store finalization
on the local executor without scratch; acquisition and scientific tasks remain
on Slurm. Actual BRAKER3 software checks passed for the selected routes. Training,
prediction output validation and biological assessment remain untested. The user
confirmed that additional AUGUSTUS predictions without hints should remain off
by default; this does not disable ES prediction or the trained models used with
evidence. Next is the real RNA+protein prediction pilot using accepted stage outputs.

13 October 2026: the user supplied the real structural_annotation README and
prediction_status.json. Sfl-n_primary completed BRAKER3 ETP training/prediction
with both paired RNA readsets and 19,445,987 evidence proteins. The ready report
records 21,528 genes, 45,399 coding transcripts, 45,399 proteins and 45,399 CDS
sequences. Structural validation passed for nonempty predictions, reference
bounds, feature IDs/parents and matching protein/CDS IDs. These are unfiltered
predictions; downstream filtering, functional annotation and biological quality
assessment remain pending. Other evidence modes have software/handoff acceptance,
not real training acceptance. Additional AUGUSTUS predictions remain disabled.

13 October 2026: start downstream implementation by reconnecting existing AGAT,
BUSCO and OMArk modules for an unfiltered baseline. AGAT retains all isoforms;
BUSCO/OMArk use a deterministic longest protein per coding gene without changing
the original predictions. Per-assembly lineage/taxid and enabled analysis settings
replace global species settings. Missing requested outputs block reporting; disabled
analyses remain explicit. The new QC-only development entry reads accepted outputs,
while main.nf uses native prediction channels. Database acquisition keeps shared
storeDir semantics; user confirmed actinopterygii_odb10 and LUCA/taxa.sqlite exist.
Python contracts passed; actual AGAT contracts and real BUSCO/OMArk execution remain
pending. MAKER AED uses available transcript/protein evidence; the legacy best-hit
protein-based CDS score is retained as an explicitly named baseline for subsequent
filter implementation, not claimed as a combined-evidence MAKER score.

13 October 2026 clarification: final AED must use all available supplied evidence,
including RNA when present. Transcript assembly (potentially StringTie) and
annotation polishing (potentially PASA or an alternative) will be discussed later;
no tool or assembly strategy is selected yet. Protein-based AED remains a historical
comparison, not the final evidence policy. Baseline contracts_01 stopped before QC
at a Nextflow 23.10.1 subworkflow declaration error. Workflow-level channel bindings
now follow the working repeat/Infernal subworkflows; corrected cluster execution
remains pending.

The real baseline run published AGAT, BUSCO and representative selection but
aborted in OMArk with exit 2. OMAmer successfully searched 21,524 proteins. The
post-analysis `omark --version` probe invoked argument validation requiring a
database; version recording now reads installed package metadata, as for OMAmer.
Existing database files are retained. The baseline report and final real OMArk
acceptance remain pending recovery with native -resume.

13 October 2026: user-supplied ready real baseline report confirms AGAT, BUSCO
and OMArk completed. Selection contains 21,528 representatives from 21,528 genes
and 45,399 coding transcripts. BUSCO actinopterygii_odb10 (3,640 groups) reports
97.9% complete, 4.0% duplicated, 0.3% fragmented and 1.8% missing. OMArk's
Percomorphaceae assessment uses 16,241 conserved HOGs: 91.71% single,
4.41% duplicated and 3.88% missing. Of 21,524 assessed proteins, 96.31% are
consistent, 0.75% inconsistent, 0% contaminants and 2.94% unknown. The four
selected proteins absent from that assessment are not assumed filtered by the
pipeline; reconcile them against OMArk classifications. Correction: the earlier
answer used the OMAmer progress count as the selection count. It is 21,528.
Real execution is accepted; this is an unfiltered baseline for subsequent
all-available-evidence scoring and filter comparisons, not a final quality claim.

The user reconciled the four omitted OMArk classifications against the selected
FASTA mapping: g11457.t1 is 2 aa; g7676.t1, g7679.t1 and g7683.t1 are 4 aa.
Each is its gene's sole isoform. They remain selected and have not been removed
by pipeline filters. The assessed count is therefore distinct from the input count;
the exact OMAmer length exclusion rule is not claimed as verified. These models
need structural inspection; short protein length alone does not establish which
legacy CDS/gene-span filters apply.

13 October 2026: user requested planning and research before reconnecting scoring
and filters. See evidence-scoring-plan.md. The proposal starts from existing
miniprot/AED, filter, TE flagging, ID mapping and QC modules; audits retained ETP
StringTie assemblies before duplicating work; adds an RNA-only assembly route;
and defines combined AED/UTR/isoform handling before implementation. PASA,
EviAnn, Mikado, junction QC and reference-transfer options are researched
candidates, not selected defaults. Refinement is evaluated on retained unfiltered
candidates before final hard filtering, with rescoring and complete sequence/ID
validation. No scientific pipeline behavior changes in this planning step.

The user accepted the proposed order and specified that canonical AED must follow
MAKER's definition and protein/transcript evidence integration. This supersedes
the proposed custom combined metric: evaluate MAKER's existing-model assessment
as the reference, require agreement tests for any compatible scorer, preserve
structures/IDs during scoring, and keep legacy protein-only/coding-region metrics
separate. Gene-level filtering aggregation is our explicit policy, not silently
claimed as part of MAKER's per-transcript AED. PASA and EviAnn are accepted for
optional refinement exploration, not promoted as default pipeline stages.

The first real evidence inventory confirms two per-readset StringTie assemblies,
merged transcript coordinates and transcript FASTA in retained GeneMark-ETP/rnaseq.
Use these as candidate reusable RNA evidence rather than the complete/incomplete/
unsupported coding subsets. Per-readset/merged BAM hints are also retained.
File presence alone does not establish compatible format/provenance; protein
alignment availability and MAKER input conversion remain to be checked. Do not
rerun StringTie just because the downstream scoring route is new.

Follow-up inventory: StringTie 2.2.1 merged output is GTF despite its .gff suffix.
It and the transcript FASTA contain the same 23,767 transcript IDs; all transcript
contigs occur in the accepted reference. Both input readsets retain reverse strand
metadata. No retained Spaln files were found. ProtHint top_chains.gff deliberately
trims CDS blocks when generating CDSpart hints and cannot substitute for untrimmed
protein alignment coverage. Reuse existing RNA evidence, create a separate cached
protein alignment with the existing miniprot module, and validate MAKER scoring-only
inputs without retraining or changing BRAKER models. The evidence inventory is
complete at the listing/ID/reference-name level; scientific conversion/runtime
acceptance remains an implementation task.

The user rejected installing full MAKER for AED assessment. Use a lightweight
implementation starting from existing interval utilities and pinned MAKER source,
not the MAKER runtime. Source revision 4bf3085ae8df4680468646a63c635094a0a3dc2a:
shadow_AED::get_AED/get_SN_SP unions supplied evidence blocks, compares them with
transcript blocks, and computes 1-(SN+SP)/2. Upstream evidence selection still
requires tracing; matching the core formula alone does not establish compatibility.
Validate against source-derived fixtures or a minimal reference scorer without
requiring a full MAKER installation. This supersedes the earlier scoring-only
MAKER runtime proposal; PASA/EviAnn remain optional refinement candidates.

13 October: implemented the first lightweight evidence AED pilot, starting from
the existing miniprot command and retained ETP/StringTie transcripts. Pinned MAKER
source separates the standard nucleotide union score from eAED and upstream
clustering. The pilot uses matching contig/strand and direct transcript-span overlap;
it does not claim to reproduce MAKER's complete clustering/preprocessing. Scores,
source availability, legacy protein/CDS comparisons and checksums are explicit.
No model filters or polishing are connected yet. All four evidence modes have
local contracts; 102 Python tests passed. A tiny real-miniprot cluster validator
checks four modes, retained RNA reuse, resume and changed-model blocking. Cluster
execution, standalone StringTie fallback and the real database pilot are pending.

The user reported passing Crest AED contracts in `testing/evidence_aed_contracts_01`
(fresh, resume and changed-model blocking). They chose the existing simpler
evidence-agreement method rather than reproducing MAKER's full clustering workflow.
Report RNA, protein and combined AED-style scores separately; do not present them
as probabilities or claim complete MAKER equivalence. Added explicit `combined_aed`
while retaining `aed` as its compatibility alias. Evaluate the real distributions
before applying the inherited cutoff; RNA junction support remains a subsequent
targeted assessment. No full MAKER port is planned.

14 October: reviewed the real AED pilot (45,399 transcripts / 21,528 genes).
82.9% of combined scores were 0.4–0.5; 98.3% had full model coverage but the
median full evidence union was 9.7 times the model length. This identifies broad
alignment-chain union as a source of compression. The user deferred Portcullis
junction support until AED selection is resolved and explicitly selected
Portcullis rather than a custom junction tool for that later increment.

Implemented a comparison selector using the closest intact RNA chain and protein
chain separately, then their union. No clipping, score rescaling or averaging is
introduced. Retain original all_union columns, selected source/alignment IDs and
chain counts. Single-chain selection is optimistic and can miss complementary
fragments within one type; scores remain agreement, not accuracy probabilities.
Keep alignments and normalization unchanged and separately cached. Validate
policy-change cache reuse, compare real distributions and inspect selected chains
before enabling removal. The revised method has not yet passed cluster or real
data acceptance; it does not promise a particular distribution.

The revised real pilot completed: combined median 0.198 (previously 0.449),
with 6.9% of transcripts at 0.4–0.5 (previously 82.9%). Best-protein scores
matched legacy scores within rounding; this is not independent validation.
The user requested evaluation of alignment/scoring usefulness before proceeding
and clarified output cleanup means removing excess AED variants, not support
fields. Keep RNA/protein/combined AED plus coverage, overlap and evidence IDs in
the main table; move legacy/broad-union comparisons and coordinate blocks to a
separate evaluation file. Remove the duplicate aed alias and production policy
switch. Cleanup does not change scores or invalidate alignment/normalization.
First evaluate coding-versus-full-transcript scope, selected alignment quality
and representative supported/conflicting/unsupported structures. Portcullis
initially complements AED after mapping; optional filtering/reassembly or use
before BRAKER requires a separate evaluated experiment. Do not blend junction
quality into nucleotide AED or turn on automatic filtering at this point.

The user clarified that legacy_protein_cds_aed should remain temporarily in the
main output during validation. Remove the identical aed alias; keep three
production scores, supporting fields and this one legacy comparison. The legacy
formula is valid for its best-protein/CDS scope, but was not a combined-evidence
or complete MAKER assessment. Equality to new protein AED in this pilot does not
guarantee equality when model CDS and exon blocks differ.

The user explicitly reaffirmed MAKER as the reference for correct AED, including
evidence preparation, and authorized a full requirements specification plus
comparison with other AED tools. This supersedes accepting the best-chain pilot
as the final method. Pinned source confirms transcript-HSP scoring rather than
an automatic protein-only CDS denominator. Added docs/aed-requirements.md with
coordinate scope, alignment/conversion eligibility, profile-dependent clustering
and model association, component reporting, reference-bag/score comparisons and
real-data acceptance. InGenAnnot, GFFtk and Panoramic are candidates to inspect,
not proven compatible substitutes. No production calculation was changed by
this requirements step; current runs remain pilot assessments and filtering is
disabled. Full MAKER prediction/polishing installation is not a production
requirement; an isolated reference harness is a development dependency.

The user clarified that a good peer-reviewed existing method is preferred over
exact MAKER equivalence and selected InGenAnnot. Skip native combined AED.
Implemented a native InGenAnnot 0.0.20 benchmark entry that uses completed models
and the same RNA/protein alignments, validates checksums and unchanged structures,
and reports matched current/native RNA/protein scores. The paper profile uses
CDS-only transcript comparison and RNA discrepancy penalty 0.25; these are not
assumed fish filtering thresholds. All parameters/process resources remain in
root nextflow.config. Native cluster acceptance is pending.

The user approved UTR refinement before final AGAT/AED/QC, with the current coding
baseline retained. Final order: BRAKER, coding baseline, native RNA-supported UTR
refinement, structural/sequence validation, final AED/QC, model filtering and
post-filter QC, then functional annotation. CDS/protein identity must be preserved
by UTR-only refinement. Native clusterize/isoform_ranking/utr_refine rank workflow
is the next integration after the unchanged-model benchmark, rather than blind
longest-transcript extension or a custom UTR algorithm. Portcullis is a separate
planned assessment; no raw BRAKER model changes are made by the benchmark.


The native AED benchmark and ten selected coordinate-level cases passed on the
real flounder data. RNA penalties can flag incomplete coding ends despite full
intron agreement, so AED-based removal is explicitly opt-in. Default
filter_max_aed and filter_aed_source are null; users may select RNA or protein
AED after inspecting distributions. The existing filter retains its gene-level
best-transcript scope and rejects unavailable requested evidence. Other structural
filter defaults remain unchanged; no filters run in the UTR pilot.

Implemented the next isolated UTR pilot in utr_refinement.nf, reusing the pinned
InGenAnnot clusterize/ranking/utr_refine workflow, accepted separate BAMs, retained
RNA or existing StringTie fallback, gffread and the existing AGAT module. Native
RNA selection is preserved. Separate BRAKER codons are retained and stop-codon
bases excluded from UTR labels. CDS coordinates/phases, extracted coding/protein
sequences, IDs/parents and reference bounds must be preserved. Original BRAKER
outputs and baseline QC remain intact. Native ranking uses stable approved input
values so report-helper edits do not force reranking. Local validation tests pass;
native cluster UTR tests and real-extension review are pending. Final native AED
rescoring and main.nf promotion follow those checks.
