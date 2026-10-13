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

RNA exons and protein alignment blocks are unioned together, including full
selected blocks outside model boundaries. Repeated hits do not add weight.
Selection requires the same contig and strand and overlapping transcript spans.
Ambiguous RNA strands are excluded and counted. This reproduces the nucleotide
core on the same selected blocks, not MAKER's entire evidence clustering,
polishing, eAED or QI. Evidence used for training is not independent validation.

All four evidence combinations are supported. With a supplied source but no
matching support, AED is 1; with neither source, AED is unavailable (`.`).
Three transcript-level AED-style agreement scores are explicit columns:
`rna_aed`, `protein_aed` and `combined_aed`. Lower is better; these are not
probabilities of annotation correctness. `aed` remains a compatibility alias for
`combined_aed`. The combined score uses the union of both evidence types, not
an average of their scores. The old best-protein/CDS comparison is separate.
Gene summaries use the lowest available transcript AED; transcript scores remain
available so this aggregation can be evaluated before filtering.

Outputs are under `structural_annotation/evidence_aed/<assembly>/`:
`aed_scores.tsv`, `aed_status.json`, `README.md`, and, when applicable, protein
alignments or new RNA assemblies. The status records reference/model/evidence
checksums, source availability, accepted/excluded evidence counts and definition.
Protein alignment verifies the genome dictionary and protein source manifest;
scoring rejects changed models or normalized evidence. Nextflow owns staging,
scratch, publication and cache. Protein alignment, RNA assembly, normalization
and scoring are separate tasks; filter thresholds do not participate in them.

## First cluster check

After syncing, use an empty test directory:

```bash
python gcl_genome_annotation/tests/validate_evidence_aed_cluster.py \
  --outdir testing/evidence_aed_contracts_01
```

The fixture uses one small protein, a retained GTF and four evidence modes. It
runs actual miniprot, checks fresh/resumed results and blocks changed predictions.
It stops after these tests. It does not run the 19.4-million-protein database or
test the StringTie fallback. Failed tests can recover with the same directory
and `--resume`.

Local acceptance covers interval arithmetic against an independent nucleotide
set oracle, evidence union rather than best-hit/averaged scoring, absent versus
unsupported evidence, opposite/ambiguous strands, UTR extensions, shared exon
parents, unchanged predictions and checksum blocking. Nextflow/cluster execution
passed on Crest in `testing/evidence_aed_contracts_01`: fresh, resume and
changed-model blocking, using actual miniprot and all four evidence combinations.
Standalone StringTie fallback and full-database performance remain pending.

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
