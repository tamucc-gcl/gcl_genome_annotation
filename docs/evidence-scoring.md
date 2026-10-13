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
Per-source scores and the old best-protein/CDS comparison are separate columns.
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
is pending; local Python tests do not establish cluster acceptance.

Once the tiny cluster checks pass, run the real scoring pilot against the
existing accepted output directory, with `--aed_inputs_root` pointing to it.
Use the same work directory and `-resume` for subsequent runs. Large miniprot
alignment retains the existing 64 CPU/300 GB/96h defaults; the fixture overrides
them to 2 CPUs/8 GB/1h. Resources are parameters in the root config. Do not change
the normal sbatch entry yet; filtering, TE flags and optional refinement remain
the following increments.
