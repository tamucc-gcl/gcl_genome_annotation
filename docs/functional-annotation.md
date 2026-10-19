# Functional annotation of the selected identified models

`main.nf` now continues from final identifiers into the existing InterProScan,
eggNOG-mapper, DIAMOND, KofamScan and DeepSig searches, then the existing functional
merge and GFF annotation logic. This increment covers protein-coding models only.
Other genome feature tracks and the unified final report remain later increments.

## Inputs and searches

The source is `structural_annotation/identified/<assembly>/`: identified GFF3,
protein FASTA, accepted transcript assessment and validation/checksums. Functional
annotation requires finalization. Existing resolved per-assembly `run_*` switches
control each tool. `--run_functional_annotation false` disables the whole increment;
its default is true. No RNA or protein prediction evidence is needed to run
functional searches on genome-only predictions.

Search-only proteins remove one terminal stop and preserve final IDs. Internal
stops, empty sequences after removing a terminal stop, and non-IUPAC amino-acid
characters are recorded as ineligible rather than silently repaired. Partial
proteins are not excluded merely for missing an initial methionine. Eligibility
never removes a model or alters the identified protein/CDS FASTAs. The existing
SeqKit chunking and search commands are reused. Per-assembly eggNOG taxonomic
scope is passed explicitly. Successful searches with no matches are valid; a tool
failure is not converted into an empty successful result.

Shared sources use their existing native `storeDir` layouts:

| Tool | Stored source |
|---|---|
| InterProScan | `interproscan_database/interproscan-data-<version>/data/` |
| eggNOG | `eggnog_database/eggnog_db/`: `eggnog.db`, `eggnog_proteins.dmnd`, `eggnog.taxa.db` |
| DIAMOND | `diamond_database/<diamond_version>.dmnd`, default `swissprot.dmnd` |
| KofamScan | `kofamscan_database/kofamscan_db/`: `ko_list` and `profiles/` |
| DeepSig | Existing tool environment and bundled models |

Absent standard sources use the existing download modules; a custom DIAMOND
index must already exist. Shared parents are created on the driver before scratch
jobs. Nextflow owns transfers and reuse. Existing databases are trusted as the
operator's declared releases; presence is not a new release/version audit. The
unversioned eggNOG/Kofam directory names retain the original layout. Do not replace
shared contents during an active run. Database inputs use symlink staging rather
than being copied into every chunk's scratch directory.

## Merge and GFF decoration

Existing parsers and product hierarchy are retained in directly testable helpers
called by `MERGE_FUNCTIONAL_ANNOTATIONS` and `ANNOTATE_GFF3`. The duplicated
structural/splice calculations, legacy AED assumptions and composite quality
logic are replaced by a direct join of the accepted assessment. RNA AED, protein
AED, structural metrics, junction support and TE flags remain unchanged.

The table has stable columns and one row per retained transcript, including
ineligible, disabled and no-hit cases. Unknown IDs and changed sources block the
handoff. Selected OMArk classes apply only to assessed representatives; other
isoforms are not assigned their representative's class.

Products retain the original priority: DIAMOND description, eggNOG preferred
name/description, significant KofamScan KO, then `hypothetical protein`. A KO below
the native threshold does not assign a product. Among multiple significant KO
hits the table selects the lowest E-value, then highest score, then KO ID. Raw
results preserve all matches. DeepSig `Chain` alone is not a signal-peptide hit.
GO terms, database references and signal/pathway notes are added without changing
coding geometry, phase, feature source or IDs. Free-text GFF delimiters are escaped.
Final geometry/ID coverage is independently compared with the identified GFF.
No further filtering or confidence score is added.

Native-format references: [InterProScan inputs](https://interproscan-docs.readthedocs.io/en/v5/InputFormat.html),
[KofamScan](https://github.com/takaram/kofam_scan),
[DeepSig output](https://github.com/BolognaBiocomp/deepsig/blob/master/deepsig/helpers.py).

## Outputs

`functional_annotation/<assembly>/` contains:

- `<assembly>.functional_annotations.tsv`: functional results and accepted assessment.
- `<assembly>.annotated.gff3`: functionally decorated coding models.
- `<assembly>.gene_product_summary.tsv`: product per gene, using the first named
  isoform in deterministic input order as the gene-level summary.
- `eligibility.tsv`, `functional_merge_status.json`, `functional_status.json`, `README.md`.
- Native merged results for each executed tool.

The Markdown README links the original identified proteins/CDS and complete maps,
structural QC and native functional results. No custom HTML is generated. Existing
structural reports remain records of the pre-functional stage; the later unified
report will connect both stages.

## Validation and execution

Local contracts cover all five formats, partial/internal-stop eligibility,
no-hit/disabled distinctions, assessment preservation, representative-only OMArk,
ID/source mismatches, product escaping and unchanged coding geometry. Nextflow
23.10.1 scripts are parsed locally. Actual scheduling and the complete five-tool
environment/database combination remain cluster/pilot checks.

First run the small cluster acceptance:

```bash
python gcl_genome_annotation/tests/validate_functional_cluster.py \
  --outdir testing/functional_contracts_01
```

It runs actual SeqKit and DIAMOND on tiny local sequences: four accepted prediction
evidence modes, an all-ineligible annotation, shared stored-index reuse by a fresh
project, no hits, disabled searches, unchanged resume and report-only cache reuse.
Native-format parser tests cover the other four tools; the fixture does not claim
to run their large databases. It runs no prediction, mapping, repeats or general
database downloads. Recovery uses the same directory and `--resume`.

After acceptance, use the normal main sbatch with `-resume` and existing database
arguments. No new config file or user checkpoint is required. A development-only
`functional_annotation.nf` entry reads completed identified/preflight/selected-OMArk
outputs with `--functional_inputs_root` for a functional-only pilot. Main passes
native task outputs directly. Search-setting changes rerun the affected tool and
merge; report-helper changes rerun reporting only. Changes to the selected set or
final numbering can legitimately require searches because query IDs/sequences changed.
