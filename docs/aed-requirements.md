# MAKER-matched AED requirements

**Historical record — not current operating instructions.** Retained for the
decision trail and flagged for archival in [cleanup](cleanup.md). Statements
below about pending stages, proposed MAKER compatibility, old AED/filter behavior
or pilot totals describe their original dates. The current implementation is
documented in [the guide index](index.md), [native AED](evidence-scoring.md),
[selection](model-selection.md) and [validation scope](validation.md). In a
conflict, use those current guides and root configuration; do not run historical
plans as additional pipeline steps.

14 October 2026. The user's target is MAKER's use of AED, including how evidence
is prepared and selected. This supersedes treating the custom best-chain pilot
as the final scoring method. This document specifies requirements, not completed
MAKER compatibility. Keep current predictions and pilot results; keep automatic
AED filtering disabled.

Subsequent user clarification: exact MAKER equivalence is not required if a good
peer-reviewed tool can be used directly. InGenAnnot is now the selected candidate;
native RNA/protein scores are benchmarked against the current method and combined
AED is deferred. See [the active InGenAnnot plan](ingenannot.md). The specification
below remains a reference/audit document, not authorization to continue a MAKER
port as the primary approach.

## Reference and equivalence target

Use MAKER source revision `4bf3085ae8df4680468646a63c635094a0a3dc2a`
as the initial fixed reference. Record its effective settings, not just a version
label. The primary route is unchanged external models (`model_gff`) plus
prealigned RNA (`est_gff`) and proteins (`protein_gff`). No gene prediction,
training or annotation polishing is required to assess external models.

There are two different equivalence claims:

1. **Scoring equivalence:** given identical models, prealigned evidence and
   effective settings, reproduce MAKER's evidence bags, coordinate blocks and
   AED values. This is the implementation acceptance target.
2. **End-to-end evidence equivalence:** using miniprot and StringTie does not
   reproduce MAKER's native alignment/polishing route automatically. Validate
   the scientific suitability and conversion of these evidence sources separately.
   Do not claim identical values across different aligners or different evidence.

Production must not require installing the full MAKER prediction pipeline.
Prefer reuse of an existing scoring implementation or isolated MAKER routines
where practical. A development-only reference harness can have dependencies
that are not required by normal pipeline execution. Decide reuse versus a port
after tracing dependencies, not by assuming a new implementation is necessary.

## 1. Model and coordinate scope

Confirmed reference behavior: `load_transcript_stats` passes the transcript hit
to `shadow_AED::get_AED`; `get_SN_SP` fills the model vector from its HSP blocks.
The standard measure compares the transcript exon footprint, including UTRs
present in that model. It does not automatically switch to CDS for protein-only
evidence. The translation offsets used for CDS/UTR reporting do not themselves
change the vector in this routine.

Let T be that model footprint, R the prepared RNA evidence blocks and P the
prepared protein evidence blocks. Following a shared prepared evidence bag:

| Production column | Model footprint | Evidence footprint |
|---|---|---|
| rna_aed | T | Union of RNA blocks retained in the bag |
| protein_aed | T | Union of protein blocks retained in the bag |
| combined_aed | T | Union of retained RNA and protein blocks |

The separate component columns are our reporting extension: restrict the
prepared bag by evidence type, rather than inventing independently optimized
component inputs. Combined AED must reproduce the reference combined score.

Proteins have no UTRs. Consequently, an excellent coding alignment can still
leave model UTRs unsupported in protein-only transcript AED. RNA can supply that
support in combined AED. A CDS-only comparison answers a different question;
retain `legacy_protein_cds_aed` temporarily and clearly label it as such.

Requirements:

- Verify whether accepted BRAKER exon and CDS footprints are identical; do not
  infer this from equal scores. Report exon/CDS lengths and UTR-capable scope.
- Preserve supplied models, feature IDs, explicit parents and alternative isoforms.
  Treat shared Parent features separately for each transcript.
- Normalize 1-based inclusive GFF and 0-based half-open internal coordinates
  exactly; test both strands and one-base boundaries. Reproduce reference handling
  of split codons, terminal stop features and model/evidence bounds.
- Distinguish absent UTR annotation from demonstrated absence of UTR sequence.
  Do not silently add UTRs, trim assembled RNA to CDS, or clip evidence to model
  boundaries to improve scores. Any reference-defined evidence transformation
  must be identified, recorded and reproduced.
- Full RNA terminal extensions may represent true UTRs or bad/fused assemblies.
  Inspect representative cases; AED alone does not determine the explanation.

## 2. Evidence alignment and eligibility

Keep the existing reference, protein databases, StringTie results and miniprot
alignments as candidate inputs. Validate identifiers, reference identity,
coordinates and source mapping before scoring. Use full RNA/protein alignment
structures, not sparse ProtHint training hints or predicted GeneMark models.

Requirements:

- Record actual aligner/assembler versions and commands, evidence sources and
  checksums. Pin accepted versions for reproducibility once audited.
- Preserve alignment metadata needed for evaluation: amino-acid identity,
  positivity/alignment score where available, aligned protein coordinates and
  query length, alignment rank, frameshifts and premature stops. Preserve RNA
  readset/library origin and strand information where available.
- Define conversion to MAKER-recognized evidence classes and HSP blocks.
  Demonstrate that miniprot GFF/PAF and StringTie GTF adapters produce the same
  genomic evidence structure intended by `protein_gff` and `est_gff`.
- Reproduce eligibility filters from the chosen MAKER input route and profile.
  Native raw-sequence alignment filters must not be assumed to apply identically
  to imported GFF evidence. Additional identity/coverage filters would change
  the method: evaluate them separately, not as undocumented compatibility fixes.
- Treat partial alignments, multimapping, homologous loci, weak alignments,
  frameshifts/stops and ambiguous RNA strands explicitly. Do not silently equate
  a genomic overlap with valid biological evidence or discard all partial
  alignments by an arbitrary cutoff.
- Missing requested input, malformed evidence or tool failure blocks assessment.
  An absent evidence type is unavailable. A supplied type with no retained local
  support has AED 1; preserve rejected-evidence reasons so this is explainable.

## 3. Evidence preparation and association with models

The existing broad-span union and closest-chain pilot are comparison methods,
not the reference selection workflow. Matching the reference requires tracing
all transformations into the final evidence bag passed to `get_AED`.

Audit and implement the applicable external-model path through `prep_hits`,
`prep_gff_data`, model/cluster association, transcript grouping and
`load_transcript_stats`. Map dependencies in `cluster`, `clean`, comparison and
hit-conversion routines. For each operation, record inputs, effective settings,
outputs, coordinate changes and included/excluded hit IDs.

The audit must cover strand separation; overlap/flank semantics; low-complexity
and redundant-splice cleaning; evidence depth handling; cluster broadening and
joining; ORF-based joining; fusion correction/splitting; single-exon, short and
unusual-splice RNA handling; and association of external models with one or
multiple evidence groups. Which operations run depends on the reference profile.
Do not assume every possible MAKER option must be implemented to match one
declared profile, but do not skip an operation that affects that profile's bag.

Profile requirements:

- Resolve and record effective values/branches for `single_exon`, `single_length`,
  `pred_flank`, `est_forward`, `correct_est_fusion` and organism type, plus internal
  clustering/splice/depth constants. Values must come from the pinned source and
  reference invocation; this document does not invent defaults.
- Preserve legitimate alternative/partial structures according to the reference,
  without creating a composite transcript by an ad hoc union of unrelated chains.
- Test adjacent, overlapping and fused loci, including contig/chunk boundaries.
  Native model IDs may differ in reference output; compare mapped structures and
  selected source IDs, not just displayed names.
- Do not rank evidence solely by whichever AED makes our annotations look best
  unless that operation is actually in the selected reference path.
- Restrict the final prepared bag by RNA/protein type for component reporting.
  No source weighting, arithmetic averaging, minimum-component shortcut, score
  rescaling or supplementary junction penalty is permitted in standard AED.

## 4. Arithmetic, reporting and acceptance

With evidence union E and overlap O = T intersect E:

`SN = |O| / |E|; SP = |O| / |T|; AED = 1 - (SN + SP) / 2`

Count each supported base once. Respect the full retained block boundaries.
Model introns are not part of its exonic footprint. Keep standard AED separate
from MAKER's eAED and quality index; they are different calculations.

Main output retains RNA, protein and combined AED, their support fields and the
temporary legacy CDS comparison. No identical `aed` alias. Detailed bag IDs,
transformations and reference comparisons belong in evaluation diagnostics.
Report source absence separately from lack of local support. Genome-only mode
has unavailable scores rather than fabricated evidence. Gene-level minimum AED
is our aggregation policy, distinct from transcript scoring.

### Acceptance gates

1. **Reference fixtures:** feed identical prealigned inputs and settings to the
   reference harness and candidate implementation. Compare accepted/rejected IDs,
   transformed HSP blocks, final bags and lengths exactly, and scores to numerical
   tolerance (before output rounding). Matching only the formula is insufficient.
2. **Structural cases:** both strands; CDS-only and UTR-bearing models; RNA-only,
   protein-only, combined and absent evidence; duplicates; competing isoforms;
   complementary partial evidence; neighbouring genes; fusion chains; unusual
   splice sites; frameshifts/stops; no local support; one-base/terminal differences.
3. **Real-data audit:** representative random and stratified flounder models,
   including g1.t1-like RNA extensions, perfect protein-only scores, source
   disagreement, short models and unsupported genes. Inspect actual sequences,
   protein alignment metadata, genomic blocks and BAM support, not just CDFs.
4. **Discrimination:** curated structures and deliberate boundary/exon/fusion
   perturbations. Explain exceptions; do not require every incorrect structure
   to have high AED when legitimate evidence for another isoform exists. Evidence
   used for BRAKER training remains evidence agreement, not independent truth.
5. **Execution:** unchanged resume, profile-only rescoring/preparation, added
   evidence/assemblies and all four evidence modes. Reuse scientific alignments
   unless their inputs/commands change. Keep root nextflow.config as the single
   configuration source and Nextflow in charge of staging/cache/publication.

Do not claim MAKER-matched AED or enable its inherited cutoff until gates 1–3
pass. Document remaining limits even after those gates. A particular histogram,
median or percentage of zero scores is not an acceptance criterion.

## Other AED implementations to compare

| Tool | Documented behavior | Relevance |
|---|---|---|
| InGenAnnot | Separate transcript/protein best-evidence scores, selectable strand requirements and transcript CDS-only scope; optional overflow penalties | Existing evaluator worth benchmarking, but its documented choices are not a claim of MAKER equivalence |
| GFFtk consensus | Pairwise transcript-coordinate AED, locus evidence association and intron-compatibility functions | Useful reusable comparison primitives; consensus selection is not automatically MAKER's bag preparation |
| Panoramic add_AED_to_gff3.py | Adds AED to EVM models using GFF evidence | Small implementation to inspect for block/input conventions; equivalence remains unverified |

Inspect pinned implementations and fixtures before adopting any as the scoring
engine. A tool using the name AED can define different model scope, evidence
selection or penalties. MAKER remains the deciding reference.

## Implementation order

1. Complete the pinned external-model caller/conversion/profile audit and freeze
   the reference harness/fixtures. Audit current model CDS/exon scope in parallel.
2. Benchmark the existing AED tools on the same small fixtures; choose reuse or
   port based on demonstrated behavior and dependency burden.
3. Implement only the evidence preparation/selection and conversion needed to
   match the declared MAKER profile, preserving existing scientific stages.
4. Run the reference, real-data and cache acceptance gates; only then replace the
   pilot as the production score and evaluate filtering thresholds.
5. Add Portcullis as complementary junction assessment. Optional changes to RNA
   assembly or BRAKER inputs are separate evaluated experiments, not hidden AED
   formula changes or prerequisites for the initial reference comparison.

## Sources

- [Pinned MAKER nucleotide scorer](https://github.com/Yandell-Lab/maker/blob/4bf3085ae8df4680468646a63c635094a0a3dc2a/lib/shadow_AED.pm)
- [Pinned MAKER preparation/caller](https://github.com/Yandell-Lab/maker/blob/4bf3085ae8df4680468646a63c635094a0a3dc2a/lib/maker/auto_annotator.pm)
- [MAKER2 AED and external-model assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC3280279/)
- [InGenAnnot AED command](https://bioger.pages-forge.inrae.fr/ingenannot/tools/AED.html)
- [GFFtk consensus API](https://gfftk.readthedocs.io/en/latest/API/consensus.html)
- [Panoramic AED script](https://github.com/MayroseLab/Panoramic/blob/master/EVM_annotation/add_AED_to_gff3.py)
