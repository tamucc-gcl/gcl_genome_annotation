# Reference and repeat preparation

`main.nf` connects normalized references, Dfam/de novo repeat preparation,
parallel masking and integrity validation to alignment/prediction. The focused
`validation/entrypoints/reference_repeats.nf` entry stops after the reference/repeat report. Both reuse
the original scientific modules with corrected contracts and native caching.

## Reference contract

Per-assembly `min_contig_length` defaults to 1,000 bp; zero disables exclusion.
Input DNA/IUPAC FASTA, including gzip, is validated for nonempty/unique records,
invalid bases, gaps and ID collisions. Unsafe first-token characters are
normalized with a reversible full-header map. Sequences are uppercased before
masking. The output contains FASTA/index, ordered contig dictionary with hashes,
exclusions and statistics. Soft masking may change case only: names, order,
lengths and uppercase sequence hashes must remain identical.

## Libraries and discovery

Default Dfam input is the complete 3.9 FamDB, all 17
`dfam39_full.<partition>.h5` files under `dfam_database/dfam_db/`.
Taxonomy resolves a supported ancestry group (for ray-finned fish,
Actinopterygii); an exact species taxid need not exist in Dfam. Explicit
`taxa_name` overrides remain authoritative and must produce nonempty extraction.
There is no runtime fallback through progressively broader groups.

Derived subsets use native stores under
`dfam_database/subsets/dfam39/tetools-1.95/taxon_<hash>/`. The subset is
checksummed; release, partition inventory/sizes and FamDB info are retained.
This is not a full checksum audit of all large HDF5 partitions. Treat stores
as immutable. RepeatPeps is reused from `dfam_database/RepeatPeps.lib`, with a
content-derived DIAMOND index under `indices/diamond-2.1.11/`.

RepeatModeler retains consensus FASTA, Stockholm, training directories and logs,
using seed 42 by default. Curation removes short, simple/low-complexity/satellite
and duplicate consensuses, with optional CD-HIT clustering. RepeatClassifier and
optional Terrier classify families. `run_terrier.py` is an active production
wrapper, not an obsolete debugging script: it applies a source/version-checked
in-memory correction for the pinned hierarchicalsoftmax single-child bug. No
installed package is edited. Preserve its pinned environment and regression checks.

The repeat-library homology filter retains RepeatPeps-supported families;
sufficient general-protein-only matches remove likely gene-derived sequences;
no-hit families remain. Main requires an existing
`diamond_database/swissprot.dmnd` or explicit `repeat_filter_database` index.
Both searches, excluded sequences and every decision remain available. Tool
failures are not converted into no-hit results. `repeat_library` supplies a
complete masking library and bypasses discovery; `repeat_reclassify=false`
disables Terrier without disabling the other repeat preparation.

## Masking and outputs

Whole contigs are balanced across at most `n_rm_chunks` batches; each occurs
once, without artificial tile boundaries. More requested chunks than contigs
does not create empty tasks. A large chromosome remains indivisible. Temporary
short-ID maps are restored safely for long names and coordinate-like suffixes.
Merged FASTA follows reference order; fragment linkage is retained and feature
IDs are unique across batches. Canonical repeat coverage uses interval unions;
per-class unions may overlap. Native `.tbl` fragment sums are not unique genome
coverage and must not be added to derive it.

| Directory | Outputs |
|---|---|
| `filtered_assemblies/<assembly>/` | FASTA, index, header map, exclusions, contig identity and stats |
| `repeats/denovo/` | Native discovery/classification and library filtering audits |
| `repeats/libraries/`, `batches/`, `logs/`, `maps/` | Combined libraries, batch accounting and native logs/maps |
| `repeats/` | Masked FASTA, GFF/OUT, alignments, summaries and native landscape |
| `reference_repeats/` | Stage status, source/reference identities, union coverage, Markdown |

Changing library filtering reuses discovery but reruns filtering/masking and
scientific descendants. Reporting-only edits need not redo repeat discovery.
Main allows RNA mapping to proceed once its unmasked reference/reads are ready;
BRAKER waits for masking and all requested evidence.

```bash
python gcl_genome_annotation/tests/validate_reference_repeats_cluster.py \
  --outdir testing/reference_repeat_contracts_01
```

Acceptance exercises actual custom-library masking, unsplit comparison, resume,
added assemblies and invalid-input blocking. It does not exercise full de novo
discovery. The real-data discovery/classification/filtering/masking run completed
and was reviewed separately. [Output guide](outputs.md), [validation scope](validation.md).
