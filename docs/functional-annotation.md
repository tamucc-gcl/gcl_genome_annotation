# Functional annotation

After final identifiers, main runs the original five scientific search tools,
merges their results and decorates the selected coding GFF3. Release packaging
then combines coding models with independent tracks. `run_functional_annotation`
defaults to true and requires finalization; per-assembly `run_interproscan`,
`run_eggnog`, `run_diamond`, `run_kofamscan`, `run_deepsig` control tools.
Genome-only predictions can still receive functional annotations.

## Eligibility, searches and stores

Search-only inputs remove one terminal stop but preserve final IDs. Empty,
internal-stop and invalid-amino-acid queries are ineligible, not silently repaired.
Missing start methionine alone does not exclude a partial protein. Eligibility
never removes a model or changes release coding FASTAs. SeqKit chunking supports
fewer actual chunks than requested and attaches settings to every chunk by key.
Native failures remain failures; successful no-hit searches remain valid.

| Tool | Existing shared source |
|---|---|
| InterProScan | `interproscan_database/interproscan-data-<version>/data/` |
| eggNOG-mapper | `eggnog_database/eggnog_db/` annotation/taxonomy DBs and DIAMOND index |
| DIAMOND | `diamond_database/<diamond_version>.dmnd`, default `swissprot.dmnd` |
| KofamScan | `kofamscan_database/kofamscan_db/ko_list` and `profiles/` |
| DeepSig | Digest-pinned BioContainers 1.2.5 image with bundled models |

Standard missing sources use native acquisition modules; custom DIAMOND indices
must already exist. Shared stores are reused with native `storeDir`, and large
database inputs use symlink staging. Do not replace store contents during runs.
The unversioned eggNOG/Kofam paths do not independently establish release identity.
DeepSig uses its configured container, not a new host TensorFlow installation;
the process initializes/checks bundled model paths and records the native version.

## Merge and decoration

The merged table has one row per retained transcript, including ineligible,
disabled and no-hit queries. Existing structural/AED/junction/TE assessment is
joined unchanged rather than recalculated. OMArk classes apply only to assessed
selected representatives. Unknown IDs, altered sources and geometry changes fail.

Product priority is DIAMOND description, eggNOG preferred name/description,
significant KofamScan KO, then `hypothetical protein`. Non-significant KO hits
do not assign products. Multiple significant KOs rank by lowest E-value, highest
score, then KO ID. Native results retain all matches. DeepSig `Chain` alone
does not count as signal-peptide support. GO/database/pathway/signal notes are
escaped for GFF3 without modifying source, coordinates, phase or IDs. Gene
product summaries use the first named isoform in deterministic input order.
There is no new composite quality score or functional-based removal stage.

## Outputs and reruns

`functional_annotation/<assembly>/` contains the functional TSV, decorated
coding GFF3, gene-product summary, eligibility, merge/status manifests, Markdown
and native merged tool results. Functional data feed the combined release and
unified report. Original identified coding sequences remain unchanged.
Changing a search setting reruns that tool and merge/descendants; report-only
changes reuse searches. A changed selected set or final ID prefix can correctly
invalidate functional queries. `validation/entrypoints/functional_annotation.nf` is a focused entry
reading completed identified/preflight/OMArk outputs via `functional_inputs_root`.

```bash
python gcl_genome_annotation/tests/validate_functional_cluster.py \
  --outdir testing/functional_contracts_01
```

Acceptance runs actual SeqKit/DIAMOND on tiny sequences, four-mode handoffs,
all-ineligible/no-hit/disabled cases, shared-store/fresh-project reuse and cache
changes. Parser tests cover all five native formats; they do not substitute for
large-database execution. All five tools completed the reviewed real-data pilot.
See [validation scope](validation.md), [databases](running.md),
[InterProScan inputs](https://interproscan-docs.readthedocs.io/en/v5/InputFormat.html)
and [KofamScan](https://github.com/takaram/kofam_scan).
