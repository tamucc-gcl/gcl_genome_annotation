# Validation record and remaining limits

Documentation status: 21 October 2026. This record distinguishes supplied cluster
acceptance/full-data reports from local tests. Historical pilot counts are not
defaults or guarantees for other assemblies. Current output manifests define the
current run; reports with `ready` establish their documented integrity scope,
not biological truth.

## Evidence available

| Area | Observed validation | Limit |
|---|---|---|
| Input/taxonomy | Four-mode/shared-RNA, malformed/blocked cluster cases and real preflight | Curated ancestry rules are finite |
| Protein preparation | Synthetic source/store/resume checks, real OrthoDB/NCBI acquisition and stored reuse | FASTA validity does not prove appropriate homology sampling |
| RNA | Actual paired/single/passthrough preparation and mapping fixtures; real paired libraries and MultiQC | Local-read pilot is not real accession/single-end validation |
| Repeats | Actual tiny masking with unsplit comparison; real discovery/classification/filtering/masking reviewed | Fixture omits large de novo discovery; no original/new numerical-equivalence claim |
| BRAKER3 | All four evidence routes in actual software-check acceptance; real ETP training and prediction | ET/EP/ES full real-data training not established by software checks |
| UTR and AED | Native tiny and connected four-mode postprocessing; real native refinement/AED and targeted evidence review | Support does not prove every UTR/coding assignment; scores are not probabilities |
| Junctions | Native low-junction fallback fixtures; real two-readset classification without fallback | Report-only support, conditional on sampled expression and alignment |
| Assessment/selection/IDs | Actual native/search/selection fixtures, cache/threshold/prefix changes and real selected set | Diagnostic flags are not automatic truth labels |
| AGAT/BUSCO/OMArk | Actual AGAT fixtures and real unfiltered/selected completeness runs | Representative proteins and native short/unclassified behavior affect denominators |
| Function | Tiny SeqKit/DIAMOND, five native parsers, shared-store/cache cases and real five-tool results | Functional matches do not independently validate model structure |
| Genome tracks | Actual native fixtures/filter switches; real filtered outputs reviewed | Telomere/centromere defaults are configurable heuristics, not functional validation |
| Release/report | Native fixture acceptance; full release checksum verification and report review | A manifest alone cannot independently reverify downloaded/remote payload bytes |
| Final sequence export | 225 local tests passed; real manifest ready with 45,060 transcripts/1,547 selected tRNAs and checksums | Separate updated cluster sequence-export acceptance summary and independent payload checksum rerun have not been supplied |

The reviewed final pilot contains 21,378 coding genes and 45,060 transcripts.
Its source before-selection set contains 21,506 genes and 45,235 transcripts.
Earlier 21,528/45,399 raw prediction counts describe an earlier run, not a
contradiction in the current before/after comparison. Do not copy those older
values into current report expectations.

## Remaining completion checks

Source consolidation passed 196 retained Python regressions, parsing of 63
Nextflow scripts/modules plus the root config, and four native-runtime
argument/routing checks. Retained interval/model helper functions were compared
with the pre-cleanup source and have identical Python ASTs. These local checks
do not replace Linux scientific-tool execution of the renamed module and
simplified preparation tuple. Use the connected postprocessing acceptance in
[the testing guide](testing.md) before the next full-data continuation.

1. Run the updated release acceptance for selected transcript/tRNA export and
   recheck exported payload checksums on the cluster.
2. Run the final chosen stage/integration regression matrix after source consolidation, including no-RNA/no-protein and disabled optional analyses.
3. Confirm ordinary `main.nf -resume` and report/selection/track-threshold changes
   preserve unaffected expensive tasks. Separate entry scopes are not promised
   to reuse all earlier standalone hashes.
4. Review source/database identity assumptions and required overrides for each
   new taxon/site. Full training/biological validation beyond the real ETP pilot
   remains a separate scope from the synthetic route matrix.

See [test commands](testing.md), [supported operation](running.md), and
[cleanup candidates](cleanup.md). No broad claim of all-organism, all-profile,
all-evidence-mode biological validation is made.
