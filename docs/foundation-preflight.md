# Input sheets and taxonomy

`main.nf` starts with the input/taxonomy preflight. `validation/entrypoints/preflight.nf` exposes the
same resolution independently and stops before acquisition or annotation.
Examples are [assembly CSV](../assets/annotation_assemblies.example.csv) and
[RNA CSV](../assets/annotation_rna.example.csv). Copy and edit them in the
consuming project; their example paths are not bundled biological inputs.

## Assembly sheet

Required columns are `assembly_id` and `assembly_fasta`; `assembly_path` is a
legacy alias. IDs must be unique and use the supported letters/digits/underscore/
dot/hyphen format. Relative file paths resolve against the sheet's directory,
not the pipeline checkout. Gzip FASTA is supported by reference preparation.

| Columns | Meaning |
|---|---|
| `taxid`, `species`, `genetic_code` | Effective organism and nuclear code; `ncbi_speciesID` is a legacy taxid alias |
| `rna_readsets` | Semicolon-separated readset IDs from the RNA sheet; empty means no RNA for this assembly |
| `use_protein_evidence` | Explicit true/false; false disables all prediction protein sources |
| `orthodb_clade`, `ref_proteome_accession`, `protein_fastas` | Taxonomy-selected OrthoDB, explicit NCBI assembly accession(s), and local protein FASTA list(s) |
| `min_contig_length` | Default 1,000 bp; zero disables length exclusion |
| `taxa_name`, `busco_lineage`, `telomere_motif`, `eggnog_taxid` | Per-assembly resolved tool settings |
| `run_busco`, `run_omark`, `run_interproscan`, `run_eggnog`, `run_diamond`, `run_kofamscan`, `run_deepsig`, `run_trnascan`, `run_infernal`, `run_telomere`, `run_centromere` | Per-assembly optional analyses |
| `individual_id`, `haplotype`, `comparison_group` | Descriptive assembly metadata; not automatic comparative analysis |

## RNA sheet

Each row has a unique `readset_id` and its `library_id`, `sample_id`, optional
`tissue`, `layout` (`paired` or `single`), plus either local `fastq_1`/`fastq_2`
or one `run_accession` (`SRR`, `ERR` or `DRR`). Do not supply both local files
and an accession. Paired local input needs both mates; single input needs one.
Readsets can be shared by assemblies without repeating preparation.

Optional row settings are `strandedness` (`auto`, `forward`, `reverse`,
`unstranded`), `trim_adapters`, `trim_quality`, `subset_fraction`, `subset_seed`.
Known protocol strand can be explicit; `auto` is resolved only after mapping.
Tables validate associations and file presence; the RNA stage validates contents.

## Precedence and taxonomy

For taxid, a nonblank row value wins over the global parameter. For species/tool
settings, explicit row values win, then explicit globals, then derivation from
the effective taxid. Row `auto` requests derivation even with a global override.
Blank inherits; `false`, `0` and `none` remain explicit. Conflicting aliases,
unknown columns and ambiguous/malformed tables are errors.

Use a fixed `taxdump_dir` containing `names.dmp`, `nodes.dmp`, optionally
`merged.dmp`. The resolver records checksums and never updates the snapshot.
Scientific name, inherited nuclear code and curated ancestry rules resolve
OrthoDB, BUSCO, Dfam and telomere settings. The rules are in
[taxonomy_rules.json](../assets/taxonomy_rules.json); OrthoDB acquisitions are
in [orthodb_sources.json](../assets/orthodb_sources.json). Mappings are finite,
not universal. Unmapped required settings need overrides or disabled analyses.
Automatic telomere coverage is initially vertebrates; no silent universal
TTAGGG fallback is applied. eggNOG search scope defaults to unrestricted `none`,
not the organism's own taxid. Close-relative proteome accessions remain explicit.

Global local-protein paths resolve against the launch directory; sheet paths
resolve against the sheet. Inspect resolved absolute paths before a large run.
The current BRAKER integration rejects nuclear genetic codes other than 1.

## Results and checks

`preflight/README.md` and `resolved_inputs.json` record readiness, evidence modes,
settings/origins, RNA associations and taxonomy/sheet/rules identities. An
assertion returns nonzero after reporting blocked requested inputs. This does
not guarantee runtime isolation from every later tool/infrastructure failure.

```bash
python gcl_genome_annotation/tests/validate_preflight_cluster.py \
  --outdir testing/preflight_contracts_01
```

Preflight cluster acceptance passed for ready four-mode/shared-readset inputs,
blocked consumers and malformed tables. Expected nonzero cases are safety tests,
not failed acceptance. Normal execution uses `main.nf`; see [running](running.md).
