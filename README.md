# gcl_genome_annotation
Pipeline for Genome Annotation

[Unfiltered annotation baseline](docs/annotation-baseline.md) reconnects the existing
AGAT, BUSCO and OMArk modules after validated predictions. `annotation_qc.nf`
assesses completed BRAKER outputs without rerunning prediction. Local contracts
pass; the real flounder baseline is ready (BUSCO 97.9% complete; OMArk 96.31%
lineage-consistent among assessed proteins). Four selected proteins absent from
OMArk's assessment were reconciled to 2–4 aa sequences retained in the input;
the synthetic acceptance summary has not been supplied.

![](assets/genome_annotation_pipeline.png)
## Pipeline rework

The new [input and taxonomy preflight](docs/foundation-preflight.md) supports two-sheet validation, global defaults, per-assembly overrides and taxonomy-derived settings. The integrated `main.nf` connects preparation, repeats, alignment, BRAKER3 prediction, optional native UTR refinement, final unfiltered QC and native RNA/protein AED. Functional annotation and model filtering have not yet been connected. Standalone entries remain available for development and testing.
Input-contract acceptance passed on Nextflow 23.10.1.

[Protein evidence preparation](docs/protein-evidence.md) adds validated,
source-preserving OrthoDB, NCBI accession and local protein inputs. Its separate
entry point stops after preparing per-assembly evidence; synthetic acceptance, real protein acquisition and same-project stored-source reuse passed.

The [rework addendum](docs/rework-decisions.md) records
end-to-end execution with native resume and existing model-filter thresholds as
the default evaluation baseline.

All stages share the root [nextflow.config](nextflow.config), which contains parameter defaults, process resources/environments, profiles and reports. Stage entry points select what runs; they do not need separate configs.

General database downloads use persistent shared stores; spike-in proteins use a project store. The real protein stored-source rerun passed; broader fresh-project reuse fixtures remain unconfirmed.

[RNA preparation](docs/rna-preparation.md) is implemented in rna_preparation.nf with paired/single FASTQ validation, raw/processed QC, conservative trimming and reproducible subsets. Synthetic cluster acceptance passed including native staging, passthrough, resume and subset-only reruns. The real local RNA pilot and separate raw/filtered MultiQC reports were confirmed by the user. It stops before alignment/repeats.

[Reference and repeat preparation](docs/reference-repeats.md) is implemented in
reference_repeats.nf, reusing the existing scientific modules with reference identity,
record-safe curation, shared database reuse and whole-contig masking repairs. Cluster
fixture acceptance passed on Nextflow 23.10.1, including actual custom-library masking,
unsplit comparison, resume, adding an assembly and invalid-input blocking. The real
discovery/classification/filtering pilot completed; its reports were reviewed as
a baseline to continue, with heavily masked small contigs flagged for later
RNA/protein support review. This entry stops before
alignment and prediction.

[RNA alignment](docs/rna-alignment.md) is implemented in rna_alignment.nf. It
uses the accepted selected reads and normalized references, produces separate
HISAT2 BAMs with mapping/strand QC, and stops before BRAKER3. There is no BAM
merge or standalone StringTie step. Python checks and actual Nextflow/HISAT2
cluster acceptance passed, including recovery, resume, strand-only changes,
adding a reference, blocked inference and no-RNA routing. The real alignment
pilot passed; both flounder readsets have validated BAMs and confidently reverse-stranded libraries. The linked alignment MultiQC report was confirmed by the user.

[BRAKER3 prediction](docs/braker-prediction.md) connects the accepted masked reference,
protein evidence and separate readset BAMs. `braker_prediction.nf` starts from
completed published stage outputs for the first real pilot, avoiding repeat discovery
and remapping. The integrated `main.nf` uses native task outputs throughout.
Input checks, runtime checks and prediction integrity reporting are implemented.
RNA+protein connected and prediction-only cluster acceptance passed on Nextflow
23.10.1, including unchanged resume and actual BRAKER3 software checks. The stage
now also routes RNA-only, protein-only and genome-only assemblies through BRAKER's
ET, EP and ES modes. Expanded four-mode cluster acceptance passed on Crest in
`testing/braker_contracts_04`, including fresh/resume runs of both entry points
and genome-only runs without optional evidence reports. The real flounder ETP
pilot passed prediction integrity validation on 13 October 2026: 21,528 genes
and 45,399 coding transcripts with matching protein/CDS records. Downstream
filtering and biological evaluation remain pending.

[Evidence AED assessment](docs/evidence-scoring.md) adds a lightweight scoring
pilot in `evidence_scoring.nf`, reusing the miniprot command and retained StringTie
transcripts. It combines supplied RNA/protein blocks using MAKER's nucleotide AED
formula, without installing MAKER or changing the predictions. Full MAKER evidence
clustering is not reproduced. Local Python acceptance passed; tiny real-miniprot
cluster acceptance and the real-data scoring pilot remain pending. Model filters,
TE flags and optional refinement are subsequent increments.

[InGenAnnot](docs/ingenannot.md) is the selected peer-reviewed evaluation candidate.
The new `ingenannot_benchmark.nf` entry compares native RNA/protein AED with the
current method using unchanged predictions and existing alignments. Combined
AED is deferred for this method. The native benchmark and selected-case review
passed on the real flounder data. `utr_refinement.nf` now provides a native
clusterize/rank/refine UTR pilot with CDS/protein preservation checks, transcript
extraction and AGAT. Tiny native UTR acceptance and the real-data UTR pilot
passed. `main.nf` now routes accepted UTR-refined models into the existing QC
modules and native InGenAnnot AED. RNA-free or refinement-disabled assemblies
retain their original models. AED filtering defaults to off; the downstream
filter module is connected after model assessment; see the current selection scope below.

The accepted model set is published under `structural_annotation/final/<assembly>/`,
with `annotation_status.json`, GFF3 and unchanged coding FASTAs. Original BRAKER
outputs remain under `structural_annotation/braker/`. Final unfiltered QC is
under `qc/final/` and native RNA/protein AED under
`structural_annotation/aed/<assembly>/`. RNA assemblies are reused from BRAKER
when available, otherwise existing StringTie assembly/merge modules supply both
UTR refinement and AED. Protein alignments are independent of UTR changes.
All downstream handoffs use native Nextflow task paths.

Before the first integrated real-data run, execute
`python tests/validate_postprocessing_cluster.py --outdir testing/postprocessing_contracts_01`
on Crest. This checks all four evidence modes with actual small UTR/AGAT/miniprot/
InGenAnnot tasks, unchanged resume, report-only cache reuse and disabled UTRs;
it performs no BRAKER training or repeat discovery. Local Python tests pass;
connected cluster acceptance remains pending. `--utr_refine false`,
`--run_baseline_qc false` and `--run_aed false` disable the respective optional
work. The default RNA AED profile assesses coding support, not UTR boundary
quality. See [the InGenAnnot guide](docs/ingenannot.md) for safeguards and scope.

[Portcullis junction QC](docs/junction-support.md) now reports exact splice-junction
support for the accepted final annotation under `qc/junctions/<assembly>/`.
Readsets are analyzed separately, with coding and UTR-associated junction summaries.
This report-only stage does not filter BAMs/models or change AED. It is enabled
by default for main runs; `--run_junction_qc false` disables it. Actual cluster
acceptance and the real-data pilot were accepted on 19 October 2026. Junction
filtering is not a mandatory pre-BRAKER step. The [continuation plan](docs/evidence-scoring-plan.md)
covers reconnecting downstream modules and adding interpretable figures across
the QC stages. `tests/validate_junctions_cluster.py` remains the small native-tool
acceptance test for future changes.

[Model assessment](docs/model-assessment.md) reconnects structural/sequence checks
and repeat-derived gene flags before selection, alongside the accepted AED and
junction metrics. It preserves all models and separates cached TE database/search
tasks from flag thresholds and reporting. Actual cluster acceptance and the real-data
pilot passed and were accepted on 19 October 2026.

[Model selection](docs/model-selection.md) follows assessment in main, adapting the
existing structural cutoffs to individual transcripts and retaining genes with a
surviving isoform. It writes separate retained/excluded GFF, protein/CDS and
assessment sets plus complete decision tables under `structural_annotation/selection/`.
AED and junction filtering remain off by default; TE flags alone do not remove models.
Original accepted models and reports remain preserved. Local validation passes;
run `python tests/validate_selection_cluster.py --outdir testing/selection_contracts_01`
on Crest before the real selection pilot. Final renaming, selected-set QC and
downstream functional annotation are not yet reconnected.
