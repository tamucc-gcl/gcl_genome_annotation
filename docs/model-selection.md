# Transcript and gene selection

`main.nf` runs selection after the accepted post-UTR assessment. This adapts the
existing `FILTER_GENE_MODELS` module and its three structural cutoffs. Original
BRAKER outputs, accepted unfiltered models, AED and junction reports are preserved.
Identifiers, coordinates, phases and coding sequences are not changed by selection.
Final renaming, selected-set QC and functional annotation are subsequent increments.

## Rules and defaults

All settings live in root `nextflow.config`. Every action accepts `off`, `flag`
or `remove`. Flags record a triggered rule without excluding the transcript.
All triggered rules are recorded, rather than stopping at the first reason.

| Rule | Default | Meaning |
|---|---|---|
| Structural | `filter_structural_action='remove'` | Apply the original three cutoffs below. |
| Single-exon CDS | `filter_min_cds_length=300` | Exclude a single-exon transcript with CDS below 300 bp. |
| Gene span | `filter_min_gene_length=200` | Exclude transcripts belonging to a gene spanning less than 200 bp. |
| All introns very short | `filter_min_intron_length=10` | Exclude a spliced transcript if every intron is below 10 bp. |
| AED | `filter_max_aed=null` | Off until a maximum and `filter_aed_source='rna'` or `'protein'` are supplied. |
| Coding junction support | `filter_min_junction_fraction=null` | Off until a minimum coding junction pass fraction is supplied. |
| TE association | `filter_te_action='flag'` | Report the assessed gene-level flag without removing models. |

Use `0` to disable an individual structural length cutoff. Equality at a minimum
length, maximum AED or minimum junction fraction retains a model. Configuring an
AED or junction threshold uses its `filter_aed_action` or `filter_junction_action`
(both default to `remove` once a threshold is supplied). These actions also allow
`flag` or `off`. No diagnostic `structure_flag` value automatically removes a model:
internal stops, partial predictions and noncanonical-splice warnings remain
diagnostic unless they independently meet one of the explicit rules above.

AED filters require an enabled AED assessment. If the requested evidence source
is absent in a particular evidence mode, that transcript's rule is explicitly
skipped rather than assigning a bad score. Junction filtering similarly requires
enabled junction QC; RNA-free, single-exon and no-coding-junction cases are
unavailable/not applicable rather than zero support. TE removal requires enabled
TE assessment. `run_model_selection=false` disables selection; otherwise main
requires `run_model_assessment=true`.

## Corrected transcript/gene semantics

The old filter summed CDS across gene isoforms and removed all children of a
flagged gene. Selection now evaluates each transcript's own CDS/exon structure and
evidence. CDS repeated in alternative isoforms is not added together. Introns
come from consecutive exons, so an absent explicit `intron` feature cannot bypass
the rule. A gene is excluded only when no transcript survives. A good alternative
can therefore rescue a gene while the failing alternative is excluded.

Explicit GFF parent relationships replace the previous current-gene and ID-name
fallbacks. Shared descendant features are written into the appropriate partition
with only surviving parents. A gene with both retained/excluded isoforms appears
in both partition GFFs, retaining its original bounds. Gene bounds can therefore
span a removed isoform; no model geometry is rewritten. Partition protein/CDS
IDs exactly match the corresponding transcript set, with no gene-name fallback.

The TE action is explicitly **gene-level**: if set to `remove`, all isoforms of a
flagged gene are excluded. This is off for removal by default because the broad
repeat-overlap/library-association flag is not definitive TE identity. The accepted
real assessment flagged 4,422 of 21,506 genes, many with strong protein and RNA
support. No new TE identity/coverage cutoff or combined confidence score is added.

## Outputs

Under `structural_annotation/selection/<assembly>/`:

- `transcript_decisions.tsv`: retained/excluded decision, every triggered rule and
  action, removal reasons, and unavailable/not-applicable enabled evidence rules.
- `gene_decisions.tsv`: transcript counts and gene decision; removal reasons list
  reasons affecting any excluded isoform, including genes retained by another.
- `retained/` and `excluded/`: synchronized GFF3, protein/CDS FASTAs and the full
  corresponding assessment columns. Original feature IDs are preserved.
- `selection_status.json`: exact settings, model/reference provenance, counts,
  reasons and output checksums; `README.md` links the results.

Empty partitions are valid; the status explicitly marks an empty retained set.
This does not establish biological accuracy or authorize downstream QC to treat
an empty annotation as a successful complete genome annotation.

## Validation and running

Local tests and Nextflow 23.10.1 parsing pass. Actual cluster runtime acceptance
and the full-data selection pilot remain pending. On Crest, first run:

```bash
python gcl_genome_annotation/tests/validate_selection_cluster.py \
  --outdir testing/selection_contracts_01
```

The test uses synthetic accepted models/assessment handoffs for all four evidence
modes and actual Nextflow selection tasks. It checks original defaults, shared
parents, isoform rescue, synchronized partitions, scientific-notation cutoffs,
source-specific AED, optional junction/TE removal, flag-only/disabled settings,
resume and threshold changes. No mapping, repeat discovery, BRAKER, native AED or
Portcullis is run. Unit tests additionally cover invalid/mutated inputs, empty
retained sets, and exact rule boundaries. Recovery uses the same directory plus
`--resume`.

After acceptance, run the normal main sbatch with `-resume`. A development-only
`model_selection.nf` entry can consume existing accepted files through
`--selection_inputs_root "${outdir}"` without launching scientific preparation.
Use the same entry and work directory for cache comparisons. Selection settings
are inputs only to the selection task, so threshold changes do not change
alignment, prediction, AED, junction or TE search task inputs.

For a completely non-removing selection preview, use:

```text
--filter_structural_action flag
```

Leave AED and junction thresholds unset and TE action at its default `flag`.
