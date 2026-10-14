# Evidence AED pilot

`evidence_scoring.nf` assesses completed predictions without rerunning BRAKER,
repeats or read alignment. It does not apply model filters or refine annotations.
This development entry is not yet connected to the normal `main.nf` run.

It starts with the existing miniprot alignment command and accepted protein
database. Current ETP runs reuse the retained StringTie merged GTF. When that
file is unavailable, or `--aed_reassemble_rna true` is requested, StringTie 2.2.1
assembles separate BAMs with resolved library strand settings and merges the
transcript annotations. It does not use the candidate GFF as an assembly guide.
The standalone StringTie route still requires cluster acceptance.

The arithmetic follows MAKER's pinned `shadow_AED::get_AED/get_SN_SP`:

`AED = 1 - (overlap / evidence_union + overlap / model_exon_union) / 2`

The first real pilot used every same-strand span-overlapping chain. Its combined
scores were compressed around 0.45: 98.3% of transcripts had full model coverage,
but the median evidence footprint was 9.7 times the model exon length. Hundreds
of differing protein chains can occupy much more space than any individual model.
This is a selection problem, not a reason to transform scores to a wider scale.

The production method selects the individual RNA chain
and individual protein chain with the lowest transcript-exon AED. Identical block
chains are deduplicated before deterministic selection; span-only matches with no
exon overlap cannot supply support. The combined score unions these two intact
chains. Full blocks remain, including extensions and incompatible exons; nothing
is clipped to the model. The formula remains unchanged.

This measures the closest available structures and is annotation-conditioned and
optimistic. It can miss complementary fragments within the same source type.
RNA and protein choices can still conflict, making combined AED worse. Exon-base
agreement does not verify junctions, reading frame or amino-acid sequence. It is
not full MAKER AED or independent biological validation. Lower or more dispersed
scores alone are not evidence that the revised metric is better.

Historical comparisons are retained in a separate `aed_diagnostics.tsv` for
evaluation, not as additional production AED columns. The production selection
switch has been removed. Best alignment/source IDs, counts, sensitivity,
specificity and overlap lengths remain in the main table. Normalization and
alignment implementations are unchanged, so `-resume` can reuse them.
Ambiguous RNA strands continue to be excluded and counted.

All four evidence combinations are supported. With a supplied source but no
matching support, AED is 1; with neither source, AED is unavailable (`.`).
Three transcript-level AED-style agreement scores are explicit columns:
`rna_aed`, `protein_aed` and `combined_aed`. Lower is better; these are not
probabilities of annotation correctness. There is no duplicate `aed` alias in
the main table. The combined score uses the union of both evidence types, not
an average of their scores. At the user's request, `legacy_protein_cds_aed` is
temporarily retained in the main table while the calculation is evaluated.
It is the best individual protein/CDS comparison, not full MAKER or combined AED.
It matches the new protein score in the current pilot but can differ when model
exon and CDS coordinates differ. The identical `aed` alias is removed.
Gene summaries use the lowest available transcript AED; transcript scores remain
available so this aggregation can be evaluated before filtering.

Outputs are under `structural_annotation/evidence_aed/<assembly>/`:
`aed_scores.tsv`, `aed_status.json`, `aed_diagnostics.tsv`, `README.md`, and, when applicable, protein
alignments or new RNA assemblies. The status records reference/model/evidence
checksums, source availability, accepted/excluded evidence counts and definition.
Protein alignment verifies the genome dictionary and protein source manifest;
scoring rejects changed models or normalized evidence. Nextflow owns staging,
scratch, publication and cache. Protein alignment, RNA assembly, normalization
and scoring are separate tasks; filter thresholds do not participate in them.

## Revised selection acceptance

After syncing, reuse the accepted tiny fixtures:

```bash
python gcl_genome_annotation/tests/validate_evidence_aed_cluster.py \
  --outdir testing/evidence_aed_contracts_01 --resume
```

The fixture uses one small protein, a retained GTF and four evidence modes. It
runs actual miniprot, checks resumed results and blocks changed predictions. It
checks that the three production AED columns and temporary legacy comparison are present alongside their support fields
and that evaluation diagnostics are published. An empty directory can instead
test fresh execution.
It stops after these tests. It does not run the 19.4-million-protein database or
test the StringTie fallback. Failed tests can recover with the same directory
and `--resume`.

Local acceptance covers interval arithmetic against an independent nucleotide
set oracle, evidence union rather than best-hit/averaged scoring, absent versus
unsupported evidence, opposite/ambiguous strands, UTR extensions, shared exon
parents, unchanged predictions and checksum blocking. Nextflow/cluster execution
passed on Crest in `testing/evidence_aed_contracts_01`: fresh, resume and
changed-model blocking, using actual miniprot and all four evidence combinations.
The original and revised full-database pilots completed and distributions were
reviewed. The revised combined median was 0.198 versus 0.449, with 6.9% rather
than 82.9% in 0.4–0.5. This verifies a distribution change, not biological
usefulness. The simplified output contract still needs cluster verification;
standalone StringTie fallback is also pending.

Once the tiny cluster checks pass, run the real scoring pilot against the
existing accepted output directory, with `--aed_inputs_root` pointing to it.
Use the same work directory and `-resume` for subsequent runs. Large miniprot
alignment retains the existing 64 CPU/300 GB/96h defaults; the fixture overrides
them to 2 CPUs/8 GB/1h. Resources are parameters in the root config. Do not change
the normal integrated pipeline entry yet; run `evidence_scoring.nf` from the
existing sbatch for this assessment. Filtering, TE flags and optional refinement
remain the following increments.

Replace the Nextflow command in the project sbatch with the following, retaining
its existing setup and `outdir`/`logdir` variables:

```bash
nextflow -log "${logdir}/.nextflow_annotation_${SLURM_JOB_ID}.log" \
  -C "${PROJECT_ROOT}/gcl_genome_annotation/nextflow.config" \
  run "${PROJECT_ROOT}/gcl_genome_annotation/evidence_scoring.nf" \
  -profile slurm -resume -w "${PROJECT_ROOT}/work" \
  --aed_inputs_root "${outdir}" \
  --outdir "${outdir}"
```

For this pilot `outdir` is `/work/birdlab/GCL/prj_seemann-flounder_genome/testing/genome_annotation`.
The run reuses the retained StringTie GTF, performs a new miniprot alignment of
the supplied protein evidence against the masked genome, computes the scores and
stops. It does not rerun repeats, read mapping or BRAKER, or remove gene models.

For the completed real pilot, rerun this same sbatch with `-resume` after the
revised tiny tests pass: only SCORE_EVIDENCE_AED should execute; alignment and
normalization should report cached. Keep the same launch directory and work
directory. The main table keeps its support fields but drops historical and
duplicate AED columns. Historical comparisons are retained in the diagnostics
file. No mathematical or evidence-selection changes are made by this cleanup.
Do not enable automatic AED removal yet.

Next accepted addition: Portcullis-based junction assessment from existing BAMs,
not a custom junction-extraction replacement. It is deliberately deferred until
this AED selection comparison is resolved.

Before that addition, follow the [alignment and scoring evaluation plan](evidence-scoring-evaluation.md).
