Questions


1. Quantify structural annotation quality (BUSCO + gene-model statistics + repeat overlap + protein-evidence support). Add OMArk for an orthology-aware “sanity check” beyond BUSCO
2. Add functional annotation and ncRNA annotation using well-supported, widely adopted tools (InterProScan + eggNOG-mapper, plus tRNAscan-SE and Infernal/Rfam).
3. Package outputs for reproducibility and submission: standardize GFF3, export FASTA/TSV/JSON deliverables, generate browser tracks, and record provenance (RO-Crate + explicit DB versions)


Functional Annotation
4. ncRNA
    - tRNAscan-SE
    - Infernal
    - Rfam

QC Functional Annotations



Other??
4. Structural “AED-like” proxies without MAKER
    - Protein support rate: fraction of predicted CDS/exons supported by protein-to-genome alignments (miniprot/Spaln3) at user-chosen thresholds (identity/coverage). miniprot supports query coverage and near-best filtering knobs to tune this.
5. JBrowse2 track
    - Assembly FASTA + index (FAI)
    - Annotation GFF3 (genes/transcripts/CDS)
    - Repeat annotation track (RepeatMasker output converted to BED/GFF3)
    - Protein-to-genome alignments (miniprot GFF3)



Repeat landscape/divergence plots
- need to output .align file
- calcDivergenceFromAlign.pl createRepeatLandscape.pl

UTR Annotation?
- can I do this with only proteins?
- 

post-prediction filtering
- no functional annotations?
- 

GFF validation/standardization

merge interpo, eggnog, diamond outputs

integrate functional annotations with gff3
- funannotate annotate

KEGG pathway mapping
- kofamscan/kofamkoala/blastkoala

signalP & TMHMM
- signal peptide/transmembrane prediction

compile qc outputs together

annotation completeness qc
- % braker proteins which got at least one functional hit

final annotation folder
- functionally annotated gff3
- protein fasta
- cds fasta
- transcript fasta
- summary stats
- tabix-indexed gff3
- samtools indexed fasta

