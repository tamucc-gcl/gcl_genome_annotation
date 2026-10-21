# Prediction protein evidence

`main.nf` validates, acquires and combines prediction protein evidence per assembly.
`protein_evidence.nf` runs the same preparation independently and stops before
RNA processing, repeats or prediction. Functional-search databases are separate
from these prediction evidence FASTAs.

## Sources

- OrthoDB odb12 partitions are selected by taxonomy/override and pinned archive
  MD5 in [orthodb_sources.json](../assets/orthodb_sources.json). An existing
  `<orthodb_database>/<clade>.fa.gz` can be checked read-only. An uncompressed
  legacy filename alone does not establish release/checksum provenance.
- NCBI assembly accessions use Datasets acquisition with `--include protein`.
  The requested assembly's exact package member must exist; the pipeline does
  not select the first arbitrary protein file. Archives, metadata and logs remain.
- Local FASTAs, including gzip, are explicit spike-ins with planning/staged
  checksums. Their original files remain preserved.

Normalization uppercases amino acids, removes one terminal stop, rejects empty
records, internal stops, invalid characters and duplicate IDs within a source,
and assigns source-specific IDs. Original identifiers/headers remain in a TSV.
Combine sources in stable source-ID order without deduplicating biological
sequences across sources. Disabled proteins yield `not_requested`, not a
sentinel masquerading as usable protein evidence.

## Stores and outputs

General OrthoDB bundles use
`<orthodb_database>/validated/<release>/normalization-v1/`. Spike-ins use
`<protein_source_store>/normalization-v1/`, project-local by default. Native
`storeDir` reuses successful sources across runs/projects using that store;
ordinary `-resume` reuses unchanged combinations. Acquisition is on workers;
store promotion uses the configured local/driver process to avoid filesystem
visibility races. Nextflow owns transport and storage placement.

Outputs include `protein_sources/<source_id>/` provenance and
`protein_evidence/assemblies/<assembly>/proteins.fa`, source mappings and
`evidence.json`; the stage README/status summarize all requested assemblies.
These are evidence sequences, not predicted gene counts. Expected acquisition
failures preserve a failed manifest/log and block consumers; unexpected tool or
infrastructure errors can terminate execution before the summary.

Preflight/source planning may rerun deliberately. Acquisition descriptors do not
include unrelated RNA/report settings. Source/network setup can be expensive;
the full vertebrate partition is not a small test. Preserve shared store contents
and use distinct directories for intended database release changes.

## Validation

```bash
python gcl_genome_annotation/tests/validate_protein_cluster.py \
  --outdir testing/protein_contracts_01
```

Tiny fixtures cover shared sources, combinations, disabled proteins, failure
isolation, adding an assembly, resume and stored-source reuse. Real OrthoDB/NCBI
acquisition and same-project reuse were confirmed. These checks establish source
integrity/transport, not biological suitability of a particular protein collection.
See [running and database paths](running.md) and [validation limits](validation.md).
