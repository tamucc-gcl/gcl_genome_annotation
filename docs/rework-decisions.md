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
