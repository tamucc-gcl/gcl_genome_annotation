# Evidence scoring, filtering and optional refinement plan

Planning snapshot: 13 October 2026. The user accepted the proposed order, specified
MAKER-compatible AED and selected PASA/EviAnn for optional refinement evaluation.
This records the agreed plan. The initial standalone nucleotide scoring pilot is
now implemented; see [implementation and acceptance scope](evidence-scoring.md).
The remaining filters, evidence assessments and refinement are planned behavior.
Start from the pre-rework modules and the accepted flounder predictions/QC.
The final AED output must use every available supplied evidence type. Transcript
assembly and refinement must support that requirement without rerunning BRAKER
solely to change scoring/filter thresholds.

## Existing components to reconnect

| Existing component | Reuse | Required adaptation |
|---|---|---|
| COMPUTE_AED / miniprot environment | Protein alignment and interval overlap calculations | Separate alignment from scoring; assembly-specific sources; add RNA; preserve legacy protein score |
| FILTER_GENE_MODELS | Structural thresholds, retained/excluded sets, reports | Missing AED unavailable; transcript-aware CDS/introns; explicit parents; multiple reasons; independent rule actions |
| FLAG_TE_GENES | Repeat overlap and translated repeat-library homology | Match accepted repeat outputs; distinguish valid no-hit output from tool failure; default flag/retain |
| RENAME_GENE_MODELS | Consistent GFF/protein/AED IDs and mapping | Extend mappings to new scores and synchronized CDS; rename after final model choice |
| MERGE_FUNCTIONAL_ANNOTATIONS / STRUCTURAL_QC | Existing sequence-integrity metrics and quality reports | Make useful structural metrics available before expensive functional searches without rewriting them |
| AGAT / BUSCO / OMArk baseline | Before/after comparisons | Distinct candidate/retained outputs and explicit representative protein sets |
| Accepted reference, separate BAMs, strand metadata, proteins and retained BRAKER working outputs | Scientific inputs | Consume native channels; validate existing ETP RNA assemblies before choosing reuse |

Keep the existing default thresholds: single-exon CDS <300 bp, gene span <200 bp,
all recorded introns <10 bp, and available AED >0.7. TE exclusion remains opt-in;
partial/internal-stop models remain flagged rather than gaining new default exclusions.
RNA-inclusive AED changes the metric: evaluate the inherited cutoff against the
protein-only baseline instead of claiming the old decisions are equivalent.

## RNA evidence route

1. Audit retained GeneMark-ETP/StringTie results: assembled per-readset transcripts,
   merged genomic exon coordinates, transcript FASTA and settings/readset provenance.
   ETP assembles each BAM with StringTie, merges transcript models, and extracts
   their sequences with gffread. Its coding projections include additional
   selection such as one isoform per locus; do not substitute the high-confidence
   training set for all supplied transcript evidence.
2. Reuse compatible retained assemblies in the current ETP pilot. Verify their
   format and reference identity; filename extensions alone do not establish GTF/GFF3.
   Audit retained protein alignments as well: full spliced alignments with source
   IDs may be reusable, while sparse BRAKER hints alone cannot reproduce the
   existing CDS-overlap protein AED. Avoid automatically launching a second
   whole-database alignment of all 19.4 million proteins without that inventory.
3. Provide a standalone StringTie route for RNA-only predictions and explicit
   RNA assembly setting changes. Use existing coordinate-sorted BAMs, one readset
   at a time, and merge transcript annotations rather than requiring a BAM merge.
   Preserve library/readset identities and support counts. Reverse-stranded
   libraries map to StringTie's --rf option; unresolved strand remains explicit.
4. Extract RNA junction and coverage support directly from BAMs as a separate
   assessment. Assembly is useful for transcript-chain comparison and refinement;
   it is not required merely to measure read coverage or intron support.
5. Compare assembled transcripts with predictions using gffcompare for exact-chain,
   partial, novel and conflicting relationships. It provides comparison labels,
   not AED or automatic annotation polishing.

Proposal: new assessment assemblies should not use the candidate BRAKER GFF as
their guide by default. A guided assembly may be useful for refinement, but its
provenance must distinguish it from an evidence assessment that did not use that
candidate annotation. RNA supporting training and scoring is evidence agreement,
not an independent validation dataset.

Sources: [GeneMark-ETP assembly and coding projection implementation](https://github.com/gatech-genemark/GeneMark-ETP/blob/main/bin/gmetp.pl),
[StringTie manual](https://github.com/gpertea/stringtie/wiki/Manual),
[gffcompare](https://github.com/gpertea/gffcompare).

### Current real-data evidence inventory

User-supplied directory listing confirms these files under
`structural_annotation/braker/Sfl-n_primary/braker_run/GeneMark-ETP/rnaseq/`:

- `stringtie/transcripts_readset01.gff` and `transcripts_readset02.gff`;
- `stringtie/transcripts_merged.gff` and `transcripts_merged.fasta`;
- `hints/bam2hints_readset01.gff`, `bam2hints_readset02.gff`,
  `bam2hints_merged.gff` and `hintsfile_merged.gff`;
- `hints/proteins.fa/complete.gtf`, `incomplete.gtf` and unsupported-class subsets.

The raw assembled per-readset/merged transcripts and FASTA are the candidate
reusable RNA evidence. Filtered coding/training classifications are supplementary,
not replacements for that evidence. BAM-derived hints can provide additional
junction corroboration, not full transcript structures.

Follow-up inventory confirms StringTie 2.2.1, two reverse-stranded input readsets
(IBC1_S123_L006 and IBC2_S124_L006), and a merged file in GTF format despite its
`.gff` extension. GTF transcript IDs and FASTA IDs match exactly: 23,767 in each.
All GTF contig IDs occur in the accepted reference dictionary. These checks
justify reusing the existing assembly rather than rerunning StringTie; full input
validation still checks bounds, strands, exon relationships and source identity.

No retained files with Spaln in their path were found. ProtHint's top_chains.gff
is unsuitable as unchanged nucleotide alignment evidence: its make_chains.py
converts CDS to CDSpart and trims exon ends by a default 15 bp, with special
handling for short blocks. Do not reconstruct full protein coverage by padding
these hints. Predicted GeneMark models and sparse hints are not protein alignments.
The proposed current-pilot route therefore reuses StringTie GTF/FASTA but creates
a separately cached protein alignment using the existing miniprot implementation,
then validates an evidence conversion for MAKER's protein_gff input. Both optional
sources feed the scoring-only assessment; no new BRAKER training is needed.

Source: [ProtHint chained-hint conversion](https://github.com/gatech-genemark/ProtHint/blob/master/bin/make_chains.py).

## Scoring contract before implementation

The canonical AED must follow MAKER's nucleotide-level calculation and its
transcript/protein evidence-handling conventions, not a newly defined composite.
The user explicitly rejected installing the complete MAKER pipeline solely for
this task. Implement a lightweight scorer using the existing interval code,
following pinned MAKER source routines and their upstream evidence selection.
MAKER is a source/reference dependency for validation, not a required runtime.
Use retained StringTie evidence and new miniprot alignments after validation;
equivalent calculation on given alignments does not guarantee equivalent values
when the alignment inputs differ.

The audited source revision is `4bf3085ae8df4680468646a63c635094a0a3dc2a`.
In lib/shadow_AED.pm, get_AED calls get_SN_SP. The core unions the provided
evidence blocks (E), unions transcript blocks (T), and measures overlap O:
SN = |O|/|E|; SP = |O|/|T|; AED = 1 - (SN + SP)/2. With no supplied hits the
routine returns AED 1. Our no-evidence-mode status remains unavailable instead
of filtering genome-only predictions on a source that was never supplied.
Audit callers/evidence selection before claiming complete MAKER compatibility;
the core routine does not itself implement all locus/strand hit selection.

Retain the legacy protein-only score separately, expose RNA agreement and junction
support, and report the canonical AED with its available evidence sources. Do not average
two AED numbers or union every overlapping isoform indiscriminately: alternative
exon chains can conflict, and transcript exons include UTRs whereas protein
alignments cover CDS. Specify locus/strand/isoform compatibility, interval scope,
UTR treatment and evidence selection must be traced to the reference implementation
before coding. Do not substitute the legacy CDS-only comparison for MAKER's
transcript-level AED. Keep any supplementary coding-region assessment distinct.
Gene-level aggregation remains an explicitly documented pipeline filtering policy,
separate from MAKER's per-transcript score.

| Supplied evidence | Scoring behavior |
|---|---|
| RNA + proteins | Both contribute; report component and combined support |
| RNA only | RNA-based support; protein component unavailable |
| Proteins only | Protein-based support; RNA component unavailable |
| Neither | AED unavailable; structural/TE checks still operate |

An absent evidence source is unavailable, not AED 1.0. Available evidence with no
matching support is a distinct condition, reported explicitly; lack of expression
in the sampled tissues is not by itself proof that a gene is erroneous. Requested
evidence acquisition/alignment failure blocks scoring rather than becoming
"unavailable" silently. A well-supported isoform should not be discarded because
another isoform lacks support; gene/transcript decisions need separate reasons.

Use the best-supported compatible transcript as the initial gene AED summary,
matching the legacy best-transcript aggregation, while retaining every transcript's
component scores. Structural rules must measure each isoform independently and
only apply the documented gene-level aggregation. Infer introns from exon chains
when needed, without mixing exons from alternative transcripts.

Fixtures must cover perfect/partial/conflicting support, both strands, alternative
isoforms, shared Parent features, absent sources, no matching support, UTR-only
differences and mismatched reference/IDs. Use pinned-source reference fixtures
and, where useful, a minimal scoring-only reference harness; the full MAKER
installation is not required. Do not claim numeric MAKER equivalence without
demonstrating it. Keep get_AED separate from the additional get_eAED routine
and QI metrics rather than inferring their behavior from similarly named scores.
Scoring must preserve every input model and ID; filtering is a separate task.

Sources: [MAKER-P evidence-specific and combined AED](https://pmc.ncbi.nlm.nih.gov/articles/PMC3912085/),
[MAKER2 pass-through assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC3280279/),
[MAKER GFF3 pass-through documentation](https://github.com/Yandell-Lab/maker),
[Pinned MAKER scoring source](https://github.com/Yandell-Lab/maker/blob/4bf3085ae8df4680468646a63c635094a0a3dc2a/lib/shadow_AED.pm).

## Refinement tools researched

These are candidates for evaluation, not promises of improved flounder annotation.

| Tool | Capability | Proposed place |
|---|---|---|
| PASA | Updates existing coding annotations using assembled transcript alignments, including exon boundaries, UTRs and isoforms | First direct refinement pilot; preserve originals and change records |
| EviAnn | Evidence-driven annotation from RNA/transcripts and proteins; accepts existing CDS annotations | Newer comparison branch; verify existing-GFF semantics rather than assuming a drop-in polisher |
| Mikado | Selects and integrates transcript/prediction models using ORF, homology and optional junction evidence | Consider when combining multiple assembly/prediction sources; additional scoring policy needs validation |
| Portcullis | Assesses/filter splice junctions from existing BAMs | Optional read-supported junction QC; BAM filtering not required by default |
| Splam | Sequence-based splice-junction scoring and annotation/alignment assessment | Modern optional junction-quality flags; validate transfer to flounder before exclusion |
| LiftOn | Transfers annotations between genomes using DNA/protein alignment | Optional comparative evidence with a suitable reference genome AND GFF; spike-in proteins alone are insufficient |

PASA remains a direct fit despite its older release history. Its official latest
release page lists v2.5.3. EviAnn has a newer peer-reviewed description in Nature
Methods (2026), but its protein-source recommendations differ from our 19.4-million
sequence BRAKER database; do not assume feeding the entire database is appropriate.
The user selected PASA and EviAnn for exploration as optional refinement routes;
neither is selected as a default yet. Additional predictors or wholesale
replacement of BRAKER are not required for this increment.

Sources: [PASA annotation updates](https://github.com/PASApipeline/PASApipeline/wiki/PASA_genome_annotation),
[PASA releases](https://github.com/PASApipeline/PASApipeline/releases),
[EviAnn inputs and behavior](https://github.com/alekseyzimin/EviAnn_release),
[EviAnn paper](https://www.nature.com/articles/s41592-026-03156-0),
[Mikado](https://github.com/EI-CoreBioinformatics/mikado),
[Portcullis](https://github.com/EI-CoreBioinformatics/portcullis),
[Splam](https://github.com/Kuanhao-Chao/splam),
[LiftOn](https://github.com/Kuanhao-Chao/LiftOn).

## Implementation and evaluation order

1. Inventory reusable ETP evidence and finalize scoring/aggregation definitions.
2. Reconnect miniprot alignment and TE flagging; normalize reusable/assembled RNA
   evidence and raw BAM junction support. Keep searches independent of filter settings.
3. Implement/validate all-evidence scoring and corrected existing structural rules;
   expose off/flag/remove actions and preserve candidate/excluded model/sequence sets.
4. Run a separately labeled refinement pilot on unfiltered candidates using PASA
   or a selected alternative. Regenerate GFF/protein/CDS/transcript sets, validate
   IDs/parents/reference/translation, and rescore changed models. Do not polish only
   hard-filtered models, since correctable candidates could already have been lost.
5. Compare raw, refined and final retained candidates with AGAT/BUSCO/OMArk plus
   RNA/protein agreement, changed exon chains, gene splits/merges, stops/partials,
   TE flags, exclusions and ID mappings. More genes/isoforms alone is not success.
6. Promote only a justified refinement policy, then reconnect functional annotation
   and final packaging. Functional tools consume the final retained model set.

Normal execution remains uninterrupted with optional analyses controlled by root
nextflow.config settings. Nextflow owns staging, scratch, store transfers and cache.
Keep separate scientific processes for transcript assembly, protein alignment,
scoring/refinement, selection and reporting. Threshold-only reruns must reuse
BRAKER, mapping, transcript assembly and protein alignments. Native -resume tests
also cover added RNA readsets/assemblies and each of the four evidence modes.

The initial planning step did not change annotation or scientific configuration.

Inventory subsequently confirmed the reusable RNA assembly and absence of raw
Spaln alignments. The lightweight scorer now has separately cached protein
alignments and existing RNA evidence. No full MAKER scoring runtime is planned.
The pinned source audit found that upstream MAKER clusters and preprocesses hits
before scoring. This pilot explicitly uses direct same-strand transcript-span
selection instead; it reproduces the nucleotide arithmetic on identical selected
blocks, not all upstream MAKER decisions. Evaluate that selection policy and the
score distributions before reconnecting AED removal thresholds.

14 October update: the real broad-union pilot compressed scores around 0.45
because selected evidence extended far beyond individual model blocks. A revised
intact-chain comparison is implemented and documented in evidence-scoring.md;
its scientific usefulness remains subject to real-data assessment. Portcullis
is the accepted next junction-support tool, deferred until selection is resolved.
Do not proceed with a full MAKER port or automatic AED filtering for this change.
