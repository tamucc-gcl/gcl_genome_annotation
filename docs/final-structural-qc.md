# Final identifiers and structural QC

After accepted model selection, main now assigns final identifiers, runs selected-set
AGAT/BUSCO/OMArk, and produces a structural report comparing the accepted unfiltered
post-UTR models with the selected identified models. Functional annotation follows as
an independent connected stage. The original accepted and excluded sets remain available.

## Final identifiers

The existing `RENAME_GENE_MODELS` scheme is retained: `<assembly>_g1`,
`<assembly>_g1.t1`, and parent-derived child IDs. Genes are numbered by reference
contig order and coordinates, with deterministic tie-breaks; transcript/child
numbering is deterministic for the same model set. Numbering can change if selection
changes. `final_id_prefix` is optional; a custom prefix is followed by the assembly
ID to avoid collisions across assemblies.

`structural_annotation/identified/<assembly>/` contains:

- `<assembly>.gff3`, `<assembly>_proteins.fa`, `<assembly>_cds.fa`: final identified
  selected models, validated for unchanged coordinates, phases and coding sequences.
- `gene_id_map.tsv` and `transcript_id_map.tsv`: every original gene/transcript,
  its selection decision/reasons and final ID. Excluded models have no final ID;
  an excluded isoform can still belong to a retained, mapped gene.
- `feature_id_map.tsv`: every retained and excluded GFF record, source partition
  and line, old/new IDs and parents, coordinates, strand and phase. Multipart
  features retain one shared final ID; anonymous retained children receive IDs.
  Shared features can occur in both partitions, with explicit parent mappings.
- `transcript_assessment.tsv`: all retained assessment columns keyed by final
  gene/transcript IDs, with original gene/transcript columns. External RNA/protein
  evidence IDs are preserved, not confused with prediction IDs.
- `annotation_status.json` and `README.md`: fingerprints, counts, source selection,
  validation and links.

Unresolved parents, malformed/mutated selections and ID collisions block publication.
An empty retained set stops finalization with an explicit selection-review error;
it is not passed to native tools as a complete empty annotation.

## Post-filter scientific QC

The existing `ANNOTATION_BASELINE` workflow and scientific AGAT, BUSCO and OMArk
modules are reused under the `SELECTED_ANNOTATION_QC` alias. Outputs go to
`qc/selected/`, leaving `qc/final/` as the accepted unfiltered comparison.

AGAT assesses all selected isoforms, including UTRs, and keeps its native statistics,
YAML and distribution figures. BUSCO runs in **protein mode**, and BUSCO/OMArk use
one longest protein per coding gene, matching the accepted before-QC policy.
For equal lengths, the original transcript-ID tie-break is preserved through the
ID map: renumbering alone cannot switch the representative protein. Each selected
representative table and exact input FASTA remain published. Native results use
the final IDs, including OMArk classification IDs.

Existing per-assembly `run_busco`, `run_omark`, resolved lineage and taxid settings
are respected. General BUSCO and OMArk databases keep their existing shared stores.
Native no-hit/short-protein behavior is retained; missing OMArk classifications are
marked unclassified on assessed representatives. Other isoforms are not assessed
by OMArk and are not assigned their representative's classification.

The original scientific helpers/commands remain unchanged for the before-QC set.
The small selected-QC wrapper changes provenance labels and preserves representative
ties. New selected model sets receive native QC; unchanged resumed tasks use
Nextflow's cache. Report-only edits do not change native tool inputs. Changing an
ID prefix changes native result identifiers and correctly invalidates selected QC.

## Structural report

Open `qc/structural/README.md` for the assembly index, or
`qc/structural/<assembly>/README.md` for the GitHub-viewable assembly report.
Pipeline-written reports are Markdown only. Existing native tool reports,
including MultiQC and RepeatMasker HTML, are retained as requested.
The report provides complete gene/transcript tables,
`structural_qc.tsv`, `completeness_metrics.tsv`, native QC summary copies and a
fingerprinted `structural_qc_status.json`.

Ten figure groups have PNG previews and exportable SVG versions:

1. Before/after gene and transcript counts.
2. Protein/CDS lengths and exon-count distributions.
3. Intron count, minimum intron length and isoform-count distributions.
4. Separate native RNA and protein AED cumulative distributions.
5. Coding and UTR-associated junction pass-fraction distributions.
6. Coding/UTR pass, fail, not-observed and unresolved-strand junction occurrences.
7. AED versus coding-junction support density plots.
8. Structural severity and repeat-association flags.
9. Explicit selection reasons and counts.
10. Before/after BUSCO completeness, OMArk completeness and consistency.

Counts, transcript/gene populations and native-tool denominators are labelled.
Junction occurrences across isoforms are not represented as distinct genomic
junctions; distinct/readset metrics remain in the linked native Portcullis report.
Missing evidence is unavailable, not zero support/AED=1. RNA and protein AED stay
separate. Retained evidence rows are re-keyed, without repeating alignment/scoring.
Gene-level repeat context is inherited from the original assessed locus and can
reflect an excluded isoform; it is not silently recomputed as new selected-only
evidence. No automatic composite high-confidence list is inferred.

The existing structural module's best/worst gene aggregation and native
BUSCO/OMArk parsers are reused/adapted. The parsers now live in `qc_metrics.py`,
separate from the retired all-in-one metric collection CLI. The new structural report uses Matplotlib
for standalone scientific figures, avoiding the legacy full functional report's
unused GO/database dependencies. Original AGAT figures are linked as well.
Native summary formats are checked before parsing: missing or malformed outputs
are not converted to zero completeness. Each QC model checksum must match its
before/after annotation. Native summaries remain available for full interpretation.

## Controls and acceptance

All defaults/resources are in root `nextflow.config`. `run_finalization=true`
enables the whole increment and requires model selection. Set it to `false` when
testing only earlier stages. Existing baseline-QC and individual BUSCO/OMArk
switches remain effective; disabling baseline QC makes the native before-comparison
unavailable, rather than inventing a baseline.

Local contracts, native before/after AGAT acceptance and the real selected-set
BUSCO/OMArk pilot passed. To check changes on the cluster:

```bash
python gcl_genome_annotation/tests/validate_final_structural_cluster.py \
  --outdir testing/final_structural_contracts_01
```

This tests actual before/after AGAT, final mappings, all four synthetic accepted
evidence-mode handoffs, ten figure groups, report source validation, resume,
report-only cache reuse, prefix changes and disabled finalization. It does not
run BRAKER, mapping, repeat discovery, AED, Portcullis or database downloads.
BUSCO/OMArk parsing and handoffs have local contracts; actual selected-set
execution was reviewed separately using the existing databases.
Recovery uses the same test directory and `--resume`.

Normal execution uses `main.nf` with `-resume`. A development-only
`final_structural_qc.nf` entry can instead read accepted published selection,
preflight and before-QC outputs via `--structural_inputs_root "${outdir}"`, without
relaunching earlier scientific stages. Keep the same entry/work directory for
cache comparisons.
