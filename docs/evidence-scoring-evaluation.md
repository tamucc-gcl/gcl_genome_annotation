# Evaluate evidence alignment and AED before filtering

The governing target is now [MAKER-matched AED](aed-requirements.md), including
prepared evidence bags. This evaluation supports that requirement; the custom
best-chain pilot does not establish MAKER compatibility.

The original broad union compressed combined AED near 0.45. Intact-chain
selection reduces that effect, but selecting the closest chain from a large
database is optimistic. A lower or wider score distribution does not establish
annotation correctness. Keep filtering disabled until evidence/model examples
show what the three scores measure.

## First assessment using existing results

1. **Check coordinate scope.** Compare model exon and CDS blocks, and selected RNA
   exon blocks. Establish whether the BRAKER models contain UTRs and whether RNA
   extensions correspond to real UTRs, adjacent genes, retained introns or fused
   transcripts. A correct coding-only model can disagree with a longer complete
   RNA transcript. The first example g1.t1 had a 3,168-bp model and 9,911-bp selected
   RNA footprint; that needs structural inspection, not automatic rejection.
2. **Check the selected protein alignments.** Use selected miniprot IDs to inspect
   alignment identity, aligned protein coverage, reported frameshifts/stops, CDS
   blocks and source provenance. The nucleotide AED alone does not encode these.
   Do not equate protein AED 0 with a correct translated protein or unique locus.
3. **Review examples across score/source patterns.** Include low and high scores,
   protein-only support, RNA-only support, RNA/protein disagreement, RNA extensions,
   alternative isoforms, short models, multi-gene chains and repeat-associated
   models. Include both random examples and extremes to avoid examining only cases
   where a chosen metric works. Inspect gene models and selected evidence together
   in genomic coordinates, and the existing BAM pileups when needed.
4. **Check discrimination.** On an explicitly separate assessment set, compare
   trusted/curated models with deliberately shifted exon boundaries, omitted exons
   and synthetic fusion models. Quantify whether scores distinguish these errors,
   not just whether they are numerically lower. Synthetic perturbations check
   sensitivity; they are not a replacement for biological validation. A model
   can be matched by an alternative alignment in a huge database.
5. **Separate agreement from validation.** The same proteins and RNA were used by
   BRAKER, so agreement is expected and not independent. If independent validation
   is needed, use curated models or evidence withheld from training; do not relabel
   the spike-in proteome as held-out if it was supplied to BRAKER.

The main score table retains three AED columns and support metrics. The separate
diagnostics table supplies selected RNA/protein blocks and model exon/CDS blocks
(0-based half-open coordinates), plus historical comparisons. These let us explain
specific disagreements without launching a new alignment or prediction run.

Required cluster files for actual alignment inspection are the accepted BRAKER
GFF3, retained merged StringTie GTF, selected miniprot GFF records and source-ID
mapping, plus the existing BAMs. Summary tables alone cannot establish identity,
frameshift correctness or the biological interpretation of terminal extensions.

## Where Portcullis contributes

Start as a complementary assessment after RNA mapping: run Portcullis on each
readset BAM and compare its observed/accepted junctions with the predicted and
assembled transcript introns. Keep library strand metadata. Running per readset
avoids Portcullis preparation merging multiple input BAMs. Do not provide the
candidate BRAKER junctions as Portcullis's reference whitelist: that could rescue
the very junctions being evaluated. Do not enable its optional BAM filtering by
default.

Report junction counts, support and quality separately from RNA/protein/combined
AED. Standard nucleotide AED does not directly measure exact intron boundaries;
blending a junction score into it would define a new metric that needs separate
validation. Single-exon models have no junction assessment, not failed support.

An optional future upstream route can use filtered alignments to reassemble RNA
evidence and then rescore unchanged models. This can change AED by changing its
evidence input; it does not change the AED formula. Filtering can lose genuine
low-expression or unusual junctions, so compare results before making it default.

Using improved RNA evidence before BRAKER is a further optional experiment and
requires rerunning prediction. Portcullis alone does not automatically alter
BRAKER BAMs or hints. Do not invalidate the accepted prediction run merely to add
junction QC.

References: [Portcullis workflow](https://portcullis.readthedocs.io/en/latest/using.html)
and [junction metrics](https://portcullis.readthedocs.io/en/latest/metrics.html).
