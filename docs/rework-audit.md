# Original-to-current pipeline audit

Source audit: 21 October 2026. The original baseline is commit `d0ebc09`
(`fix scaff name bug`, 24 August 2026), the parent of `5c68a28`
(`inital input rework`, 9 October 2026). The current baseline is `f671e90`.
This compares workflow connections, helpers, defaults and output contracts; it
does not rerun the original pipeline or establish numerical equivalence.
The supplied full-data `annotation_report.md` was also checked as evidence of
the current output layout, not as an independent validation of every source file.

Follow-up: final selected transcript and tRNA sequence export is now implemented
in release packaging, with local checks for both strands, exon splicing, selected
IDs, disabled/no-hit tRNA tracks and checksums. The inventory below records the
gap as found at the audited commit. Cluster acceptance of this addition remains
to be run; see [release packaging](annotation-release.md).

## Overall result

The major scientific stages in the original `main.nf` remain connected in the
current pipeline. Most original scientific tools are reused, while orchestration,
validation and several helpers have changed substantially. Retaining a module
name does not mean its implementation or scientific results are identical.

The main remaining output-parity gap is final sequence export: selected,
final-ID spliced transcripts and selected tRNAs are not included in the release.
Documentation and obsolete development paths also need a cleanup pass.
No production code or default was changed as part of this audit.

## Stage comparison

| Capability | Original | Current | Assessment |
|---|---|---|---|
| Inputs and taxonomy | Assembly sheet and mostly global defaults | Validated assembly/RNA sheets, resolved taxonomy and per-assembly settings | Replaced and extended |
| Reference preparation | Filter assembly | Filter/normalize reference, contig dictionary, sequence identities and checks | Retained with stronger contracts |
| Protein evidence | OrthoDB plus reference-proteome preparation | Validated OrthoDB, NCBI and local sources; shared/project stores and per-assembly combinations | Replaced preparation; retained purpose |
| RNA preparation and mapping | No RNA path in original main workflow | Acquisition/validation, raw/processed QC, trimming/subsetting, HISAT2, mapping and strand QC | Added |
| Repeat discovery | RepeatModeler, curation, RepeatClassifier, Terrier | Same scientific chain, with interface and correctness fixes | Retained and adapted |
| Repeat-library filtering | RepeatPeps and non-TE protein homology | Same two-source purpose, explicit index/provenance and validated libraries | Retained and adapted |
| Genome masking | Parallel RepeatMasker | Parallel RepeatMasker, whole-contig batches, validated merge and reference integrity | Retained and corrected |
| BRAKER | Softmasked genome and optional protein argument; no RNA input | BRAKER3 ETP/ET/EP/ES routing with separate BAMs and validated evidence | Extended; all four combinations explicit |
| UTR refinement | Absent | Optional InGenAnnot refinement with rejected-proposal audit and coding-sequence invariance checks | Added |
| AED | Custom protein-derived score and no-protein stub | Native InGenAnnot RNA and protein scores; unavailable evidence stays unavailable | Intentionally replaced, not numerically interchangeable |
| Junction evidence | Sequence-based splice checks within downstream metrics | Retained sequence checks plus separate Portcullis support assessment | Extended |
| Model assessment and selection | Structural/AED filtering, then TE flagging | Explicit transcript assessment, gene-union TE evidence, independent opt-in rules and complete decisions | Reorganized; structural thresholds retained |
| Final IDs | Rename retained coding models | Deterministic final IDs with complete gene/transcript/feature maps and excluded decisions | Retained and extended |
| Structural QC | AGAT, BUSCO, OMArk and structural summaries | AGAT/BUSCO/OMArk before and after selection, explicit representative provenance and structural figures | Retained and extended |
| Functional searches | InterProScan, eggNOG, DIAMOND, KofamScan, DeepSig | All five retained with eligible-input preparation, per-tool switches and native result merging | Reconnected and adapted |
| Functional decoration | Merge functional tables and decorate combined GFF3 | Decorate selected coding models, then package with separate genome tracks | Retained; packaging separated |
| tRNAs | tRNAscan-SE on softmasked reference | tRNAscan-SE on filtered unmasked reference, native high-confidence filtering and raw/selected audits | Retained; input and selection changed |
| Other ncRNAs | Parallel Infernal/Rfam on softmasked reference | Same search, native inclusion/overlap selection and raw/selected audits | Retained and extended |
| Telomeres | tidk windows and terminal/internal labels | Same search plus configurable per-window density and actual-end-distance selection | Retained and extended |
| Centromeric repeats | centroAnno assembly mode and GFF conversion | Same assembly-mode purpose; pinned executable build, coordinate validation and configurable quality selection | Retained and corrected |
| Coding protein/CDS FASTAs | Extract after final functional decoration | Preserve validated BRAKER coding sequences through selection and ID assignment | Delivered, but header behavior differs |
| Spliced coding-transcript FASTA | Final extraction via gffread | UTR stage publishes unfiltered `transcripts.fa`; selected final-ID transcript FASTA is absent | Output gap |
| Selected tRNA FASTA | Final extraction from merged GFF3, including exon splicing | Native raw tRNAscan FASTA remains; no export synchronized with selected tRNA GFF3 | Output gap |
| Final release | Combined annotation and sequence outputs spread across stages | Combined/separate GFF3s, references, coding FASTAs, ID maps, tables, provenance and checksums | Added coherent packaging |
| Unified report | Markdown plus R figures, confidence summaries and tuning guidance | Markdown plus Matplotlib/reused QC figures, explicit populations/settings and source manifests | Replaced presentation; results/methods-only policy intentional |

## Intentional behavior changes to document

- The default minimum retained contig length changed from 50,000 bp in the
  original `main.nf` to 1,000 bp in the root configuration. The reference entering
  every downstream stage can therefore differ even with the same input FASTA.
- Single-exon CDS <300 bp, gene span <200 bp and all recorded introns <10 bp
  remain the structural baseline. Correct transcript aggregation and explicit
  parents mean old and new decisions need not be identical.
- The original AED cutoff was 0.7; current AED selection defaults to off and
  requires an explicit RNA/protein source and cutoff. The scoring method also
  changed. Combined/legacy AED is not a production metric.
- TE evidence defaults to flagging, not removal. Junction-based removal is
  optional. Missing evidence is distinct from a score of zero or one.
- `braker_ab_initio=false` remains the default for the extra AUGUSTUS
  no-hint output. Genome-only routing still explicitly uses ES; these are
  different controls.
- tRNAscan receives the filtered unmasked reference now; Infernal retains the
  softmasked input. New native/confidence filters mean the selected tracks are
  intentionally smaller than the raw searches.
- BUSCO and OMArk assess explicit representative proteins rather than treating
  all alternative isoforms as independent loci. Before/after uses consistent
  populations; it is not directly equivalent to every old summary.
- Final protein/CDS FASTAs currently contain final identifiers alone. Products
  are delivered in the functional GFF3 and linked tables. The old extraction
  module was intended to produce descriptive sequence headers; restoring that
  convenience requires an explicit final-export decision.
- The old synthetic high-confidence label and threshold-based interpretation
  are not restored in the unified report. Structural diagnostic flags remain;
  no new composite confidence score is calculated.
- Native HTML reports remain. Pipeline-written HTML reports are excluded.
- Current release manifests define exported deliverables. Older published files
  can remain after settings change; they do not become current results merely
  because they still exist on disk.

## Concrete follow-up work

### 1. Finish final sequence export

Use selected final-ID coding GFF3 and selected tRNA GFF3 as the source of truth.
Reuse the original extraction logic where compatible, rather than reconnecting
the old process unchanged. It selects coding records by a hard-coded `BRAKER3`
source column and removes out-of-bounds models during extraction, both of which
conflict with the current explicit-model and fail-on-invalid-input contracts.
Current native tRNA annotations also require checking against the old hierarchy
assumptions before reusing its splice-extraction code.

The output should synchronize transcript IDs with the selected coding annotation
and tRNA IDs with the selected track, handle both strands/introns, include UTRs
in coding-transcript sequences, preserve coding protein/CDS content, and add the
new deliverables to release manifests/checksums. Disabled/no-hit tRNA tracks
must retain their explicit status. Descriptive FASTA headers are a separate
presentation choice and must not change identifier tokens or sequence contents.

### 2. Reconcile documentation with the accepted implementation

The root README still describes report acceptance, track packaging, evidence
scoring, downstream evaluation and some shared-store testing as pending even
though later sections and supplied acceptance runs establish completion. It
also contains historical pilot counts alongside current workflow descriptions.
Replace this chronological account with current usage, supported modes, defaults,
outputs and limitations. Apply the same distinction to development-stage docs.

### 3. Classify superseded code before retiring it

Candidates not used by the production `main.nf` include the old protein-source
modules (`downloadOrthoDB.nf`, `prepareProteinDb.nf`), custom protein AED
(`computeAed.nf`), the alternative `repeatMask.nf`, and the R report plotting
script. The original `extractSequences.nf` is unconnected but contains logic
needed for the sequence-export gap; retain it until that work is settled.

The custom `evidence_scoring.nf` and `ingenannot_benchmark.nf` remain development
entries. Their names in the repository do not mean that production still runs
the old scoring or benchmark outputs. Decide which development interfaces and
tests remain supported before removing their helpers/configuration labels.
For example, configuration still includes old extraction/protein-preparation/AED
labels. They must be checked against the chosen supported entry points, not
deleted solely because production no longer invokes them.

### 4. Perform final regression after the output gap is resolved

The supplied cluster acceptances cover individual stages, four-mode synthetic
handoffs and cache changes; the full-data run exercises RNA+protein execution.
They do not establish full real-data BRAKER training for every evidence mode or
bit-for-bit parity with the original pipeline. The final regression should cover
the supported end-to-end mode/disable matrix and ordinary resume, and confirm
that export, threshold and report changes preserve unaffected expensive results.
Check final sequence/GFF/table IDs, complete maps and release checksums together.

## Audit boundaries

This inventory finds no missing major scientific search stage. It does not
certify all algorithms, every parser field or all taxonomic groups. For example,
the current BRAKER planner explicitly rejects unsupported nonstandard nuclear
genetic codes rather than silently treating them as the standard code. Such
limits belong in the supported-scope documentation, not in claims of universal
organism support.

No new scientific tools, filters, confidence scores or output deletions are
proposed by this audit. The immediate implementation priority is restoring final
sequence-export parity, followed by documentation and the agreed regression pass.
