# Tests and development entries

Normal use is `main.nf`. Keep tests for active scientific behavior, input/output
contracts and native cache handling even when the initial rework milestone passed.
Past acceptance is evidence for that revision, not a replacement for regression
after later changes. [Validation scope](validation.md) distinguishes real tools,
synthetic accepted handoffs, software checks and real-data training.

## Local Python suite

From the pipeline checkout with Python 3.11 and Matplotlib available:

```bash
python -m unittest discover -s tests -p 'test_*.py'
```

Imports used by tests add repository helper paths through their fixtures. Plot
tests can skip without Matplotlib; read the test summary before claiming complete
coverage. Most contracts use stdlib/synthetic data; the native cluster validators
exercise Nextflow/tool behavior separately. The retained local suite passes 196
tests after source consolidation. Superseded custom
arithmetic/comparison and one-off investigation tests were removed; native AED,
model invariance, selection, release and reporting regressions remain. Test counts
can change as obsolete coverage is retired; preserve useful behavior coverage. The Python test modules themselves share fixtures, so deleting
an older-looking module can break retained tests.

## Cluster acceptance

Run controllers from the consuming project using its usual Nextflow/Conda/
Singularity setup. Choose a new empty test outdir; `--resume` recovers an
interrupted validator using its existing inputs. Stop the old controller first.
Each validator keeps `acceptance.json`, console/Nextflow logs, trace and synthetic
fixtures. Some deliberately run commands expected to fail and verify the
failure report. The validator's final status, not one nested exit code, determines
acceptance. Many validators also run the entire Python suite before native tasks.

```bash
python gcl_genome_annotation/tests/validate_release_cluster.py \
  --outdir testing/release_sequence_contracts_01
```

| Validator under `tests/` | Native scope; important exclusions |
|---|---|
| `validate_preflight_cluster.py` | Resolution, four-mode/shared-RNA and failed-input contracts; no annotation |
| `validate_protein_cluster.py` | Tiny local protein sources, stores/combination/resume; real network acquisition is separate |
| `validate_rna_cluster.py` | FastQC/fastp, paired/single/passthrough/subset/no-RNA; no mapping |
| `validate_reference_repeats_cluster.py` | Actual custom-library masking and split/unsplit comparison; no de novo discovery |
| `validate_rna_alignment_cluster.py` | Actual HISAT2/BAM/strand checks and cache; no prediction |
| `validate_braker_cluster.py` | Connected preparation/alignment and four-mode BRAKER software checks; no training |
| `validate_annotation_baseline_cluster.py` | Actual AGAT, representatives, integrity/cache; no full BUSCO/OMArk database run |
| `validate_utr_cluster.py` | Native cluster/rank/refine, gffread/AGAT, both model strands, no-RNA/disabled/cache |
| `validate_postprocessing_cluster.py` | Connected UTR/QC/miniprot/native AED for four accepted-mode handoffs; no BRAKER training |
| `validate_junctions_cluster.py` | Actual Portcullis low-junction fallback/readset/model support; no remapping |
| `validate_assessment_cluster.py` | Actual SeqKit/DIAMOND/BEDTools assessment, no-hit/disabled and cache |
| `validate_selection_cluster.py` | Rule actions, partitions, source-specific AED/junction/TE options; no scientific searches |
| `validate_final_structural_cluster.py` | IDs/maps, before/after AGAT, reporting/prefix/cache; large BUSCO/OMArk separate |
| `validate_functional_cluster.py` | SeqKit/DIAMOND, five-format parsers and fresh-project stores; four other full tools separate |
| `validate_genome_tracks_cluster.py` | tRNAscan/high-confidence filter, synthetic calibrated Infernal, tidk/centroAnno, selected/raw/cache |
| `validate_release_cluster.py` | Packaging, final transcript/tRNA FASTAs, IDs/maps, checksums, function/track switches and cache |
| `validate_report_cluster.py` | Nextflow/Matplotlib report, two releases and report-only cache; no searches |

The old custom-score/comparison validators have been removed with their obsolete
entries. Native AED remains exercised by `validate_postprocessing_cluster.py`
and the retained production score/invariance tests. Shared model/native-score
fixtures are in `annotation_fixtures.py`, rather than imported from a benchmark
test class. See [consolidation record](cleanup.md).

The four `check_*.groovy` files test production argument construction/routing
with a Nextflow runtime. Keep them alongside native validators. Every current
`*_fixture.nf` and `utr_fixture_bams.nf` supports an active stage validator.
Those tiny workflows are tests, not extra production pipelines.

## Focused root entry points

These are retained diagnostic/test interfaces, not required manual checkpoints.
Published-input readers validate source identities and require backing inputs
to exist. Switching entries or workflow scopes can change cache identities.

| Entry | Stops after / accepted-input purpose |
|---|---|
| `preflight.nf` | Sheet/taxonomy resolution |
| `protein_evidence.nf` | Protein preparation |
| `rna_preparation.nf` | RNA preparation |
| `reference_repeats.nf` | Reference/repeat preparation |
| `rna_alignment.nf` | Mapping/QC using accepted reads/reference |
| `braker_prediction.nf` | Prediction using accepted preparation/alignment |
| `annotation_qc.nf` | Raw BRAKER QC |
| `utr_refinement.nf` | UTR refinement and standalone QC |
| `junction_qc.nf` | Model junction assessment |
| `model_assessment.nf` | Assessment using completed scientific results |
| `model_selection.nf` | Selection using accepted assessment |
| `final_structural_qc.nf` | Final IDs and selected QC |
| `functional_annotation.nf` | Functional annotation using accepted identified models |

The custom scoring and benchmark root entries are retired. Moving retained
diagnostic entries under a development directory is a possible
later reorganization, not a safe rename alone: includes, `projectDir`, helper
paths, docs and validators must change together. Those retained diagnostic entries
were not relocated. Remaining proposed archives are listed separately. Keep actual project work, reports and environments separate
from pipeline source; never remove active cache solely to tidy the checkout.

## Source-consolidation regression

After syncing the cleanup, use a new directory for the native connected check:

```bash
python gcl_genome_annotation/tests/validate_postprocessing_cluster.py \
  --outdir testing/postprocessing_consolidation_01
```

This exercises the renamed native AED module and simplified evidence-preparation
tuple. Then use normal main resume. Process names and scientific commands remain
stable, but renamed staged helper files change some hashes: a one-time rerun of
affected UTR/AED/assessment/QC descendants is possible. Do not assume this cleanup
is report-only or force old cache hits for changed helper identities. Mapping,
repeat discovery and BRAKER scientific inputs were not intentionally changed.
