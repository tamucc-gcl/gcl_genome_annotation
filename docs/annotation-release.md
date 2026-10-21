# Final annotation packaging

`main.nf` packages accepted identified models, optional functional
decoration, RepeatMasker annotations and enabled **selected** genome tracks.
It consumes native task outputs, not asynchronous `publishDir` files. It performs
no searches, model selection, renaming of coding models or new scientific filtering.
The existing finalization switch controls whether identified models and their
release are produced. BRAKER software-check-only runs omit packaging.

## Deliverables

Start at `release/README.md`, then `release/<assembly>/README.md`.

| Path within each assembly | Meaning |
|---|---|
| `<assembly>.gff3` | Combined coding, repeat and enabled selected genome tracks |
| `<assembly>.coding.gff3` | Coding models alone, with functional attributes when enabled |
| `<assembly>_proteins.fa`, `<assembly>_cds.fa` | Exact identified selected coding sequences; IDs unchanged |
| `<assembly>_transcripts.fa` | Spliced selected coding transcripts, including accepted UTRs, with final transcript IDs |
| `<assembly>_tRNA.fa` | Spliced selected tRNAs with selected-track IDs; omitted when the tRNA track is disabled |
| `<assembly>.reference.fa`, `<assembly>.softmasked.fa` | Normalized retained reference, unmasked and softmasked |
| `tracks/<assembly>.<track>.gff3` | Separate repeat/tRNA/ncRNA/telomere/centromeric-repeat annotations |
| `id_maps/` | Existing complete gene/transcript/feature maps, including excluded decisions, plus a track ID map |
| `tables/` | Identified transcript assessment; functional and gene-product tables when enabled |
| `provenance/` | Accepted reference dictionary and source-stage manifests |
| `release_status.json`, `checksums.sha256` | Current contents, counts, reference identity and file checksums |

Coding IDs and sequence contents are preserved. Selected genome-track IDs are
preserved. RepeatMasker IDs gain an assembly-specific `<assembly>.repeat.` prefix
to avoid cross-track collisions, with every record represented in `track_id_map.tsv`.
External evidence identifiers, such as a repeat's `Target`, remain unchanged.
The combined file groups tracks, emitting parents before children; it is not a
globally coordinate-sorted file. It retains overlaps between evidence tracks.

Transcript and tRNA FASTAs are extracted from the accepted final annotation and
the unmasked retained reference. Exons are joined in transcript orientation on
both strands; tRNAs without exon children use their annotated span. Transcript
exports include all selected exons, including accepted UTRs, rather than only
CDS blocks. They do not include removed models or candidates rejected by the
tRNA confidence filter. Headers contain final identifiers only, with products
available in the functional tables/GFF3. Existing protein/CDS files are copied
unchanged. Extraction performs no prediction or filtering and blocks unresolved
strands or invalid/overlapping exon blocks instead of silently skipping models.

An enabled tRNA track with no selected tRNAs produces an empty tRNA FASTA;
a disabled tRNA track produces none. `sequence_counts` in the release manifest
records exported transcript and enabled-track tRNA sequence counts. Both FASTAs
are included in the normal output checksums. This export change requires only
packaging and downstream index/report tasks to rerun under normal resume.

An enabled search with no hits has an empty GFF3 track. Disabled tracks are
recorded explicitly and have no current deliverable. Functional annotation can
be disabled: packaging then uses the accepted identified coding GFF3 and omits
functional tables. There is no additional packaging configuration file.

The release is a self-contained set of scientific deliverables, copied by native
Nextflow publication so that work-directory cleanup does not break it. Detailed
QC reports, native tool reports, raw candidates and rejection audits remain in
their established stage directories, linked from the release README. No files
from previous runs are deleted automatically. Existing published files can
remain after optional stages are disabled; `release_status.json` and its checksum
list define the **current** release. Copy only listed files when exporting it.
The checksum list includes the manifest; the manifest lists all payload files
but excludes itself and the checksum list to avoid circular hashes.

## Safeguards

Packaging blocks altered identified files, mismatched references or masking,
missing enabled tracks, changed selected-track checksums, invalid coordinates,
inconsistent multipart IDs, cross-track ID collisions, dangling/cyclic parent
relationships, mismatched coding FASTA/table IDs and incomplete final ID maps.
Functional GFF3 geometry is checked using the existing functional validation
helper. No ncRNA or repeat feature is counted as a protein-coding gene.

The final release does not replace the original scientific outputs or imply
that centromeric repeat predictions establish functional centromeres. The
[unified Markdown report](unified-report.md) is a separate downstream stage.

## Acceptance

After syncing, run on Crest with the usual Nextflow/Conda modules:

```bash
python gcl_genome_annotation/tests/validate_release_cluster.py \
  --outdir testing/release_contracts_01
```

This exercises actual Nextflow packaging for two accepted synthetic assemblies,
functional annotation on/off, enabled empty and disabled tracks, synchronized
coding outputs/maps/checksums, resume and index-only changes. It performs no
prediction, mapping, repeat discovery or functional searches. Use `--resume`
with the same directory after a failure.

Then resume the usual `main.nf` run. New packaging tasks should run while existing
scientific tasks remain cached. Review `release/README.md` and each assembly's
`release_status.json` before exporting the release. Cluster packaging acceptance
and the original full-genome release review passed; all original payload checksums
were verified. The subsequent real sequence-export manifest is ready with 45,060
transcripts and 1,547 selected tRNAs. The updated export has 225 passing local
tests; its separate cluster acceptance summary/checksum rerun has not been
supplied. See [validation record](validation.md).
