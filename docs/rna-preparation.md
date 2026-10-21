# RNA preparation and QC

`main.nf` prepares each requested readset once before reference-specific mapping.
`validation/entrypoints/rna_preparation.nf` is the focused preparation entry. Local reads remain declared
path inputs; accession inputs are acquired with SRA Toolkit. Nextflow manages
staging/scratch, without helper-created passthrough links or copies.

## Content validation and defaults

The stage checks four-line FASTQ records, gzip readability, nonempty sequence/
quality, base/quality validity, read IDs, paired counts/order and mate identity.
It records SHA-256 identities and read units: pairs for paired data, individual
reads for single data. Wrapped FASTQ is outside the supported format; quality
range checks alone cannot distinguish every legacy encoding.

| Control | Default behavior |
|---|---|
| `trim_adapters` | true; fastp paired-end adapter auto-detection |
| `trim_quality` | true; 3-prime tail trimming, four-base Q20 window |
| `rna_min_read_length` | 30 bases after trimming |
| `rna_max_n_bases` | 5 |
| `subset_fraction`, `subset_seed` | 1.0, 42 |

Both trimming flags false means pass through existing reads without fastp
filtering. When trimming runs, general low-quality-percentage filtering and
automatic polyG trimming are disabled. No deduplication, read merging, orphan
rescue or polyA trimming is introduced. Subsetting follows full cleaned-read QC;
SHA-256 decisions on seed/read ID keep mates together. The fraction is expected,
not an exact count. Fraction 1 reuses paths; an empty subset blocks consumers.

FastQC 0.12.1 assesses full raw and processed reads; fastp 0.24.0 keeps JSON/HTML
and retention logs; MultiQC 1.25.1 separates raw and filtered reports. FastQC
warnings are reported, not automatic readset exclusions. `auto` strand remains
pending until reference-specific alignment; preparation-ready does not resolve it.

## Outputs

| Under `rna/` | Contents |
|---|---|
| `raw/<readset>/` | Validation/source manifest; downloaded reads for accessions |
| `processed/<readset>/` | Manifest, new cleaned reads and fastp results when enabled |
| `selected/<readset>/` | Manifest and new subset reads when requested |
| `qc/raw/`, `qc/processed/` | Native FastQC outputs |
| `multiqc/raw/`, `multiqc/filtered/` | Separate native reports and data/logs |
| `README.md`, `rna_status.json` | Per-readset and per-assembly readiness, actual selected input paths |

Native reports remain HTML where produced by the tools. No-RNA runs identify
the absence of inputs. Published new read bundles may link to work: retain work
and original FASTQs. Expected content/acquisition failures are reported and block
consumers; QC/tool/infrastructure failures remain task failures. There is no
general shared RNA store. Subset-only changes reuse acquisition, trimming and
full-read QC; reporting edits do not reprocess reads.

```bash
python gcl_genome_annotation/tests/validate_rna_cluster.py \
  --outdir testing/rna_contracts_01
```

Actual-tool paired/single/shared-readset, passthrough, subset-change, failure,
no-RNA and resume acceptance passed. The real local paired-end preparation and
separate MultiQC reports were accepted. Synthetic single-end support does not
constitute a real single-end biological pilot; real accession-download coverage
is distinct from the validated local-read pilot. See [mapping](rna-alignment.md).
