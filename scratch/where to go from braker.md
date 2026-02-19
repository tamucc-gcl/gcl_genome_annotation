Questions


1. Quantify structural annotation quality (BUSCO + gene-model statistics + repeat overlap + protein-evidence support). Add OMArk for an orthology-aware “sanity check” beyond BUSCO
2. Add functional annotation and ncRNA annotation using well-supported, widely adopted tools (InterProScan + eggNOG-mapper, plus tRNAscan-SE and Infernal/Rfam).
3. Package outputs for reproducibility and submission: standardize GFF3, export FASTA/TSV/JSON deliverables, generate browser tracks, and record provenance (RO-Crate + explicit DB versions)



QC Structure
1. Busco on proteins
2. OMArk - gene repertoire consistency/contamination
3. Gene model descriptive stats - AGAT 
    - Gene length, CDS length, exon count per transcript, exon length, intron length
    - Fraction of single-exon genes
    - Start/stop codon completeness (if tracked)
    - Repeat overlap statistics (CDS and gene span)
    - `agat_sp_statistics.pl --gff annotation.gff3 --output agat_stats`

Functional Annotation
1. InterProScan
    - genome scale functional classification
2. eggNOG-mapper
    - go/kegg terms
3. diamond
    - gene names
4. ncRNA
    - tRNAscan-SE
    - Infernal
    - Rfam



Other??
4. Structural “AED-like” proxies without MAKER
    - Protein support rate: fraction of predicted CDS/exons supported by protein-to-genome alignments (miniprot/Spaln3) at user-chosen thresholds (identity/coverage). miniprot supports query coverage and near-best filtering knobs to tune this.
5. JBrowse2 track
    - Assembly FASTA + index (FAI)
    - Annotation GFF3 (genes/transcripts/CDS)
    - Repeat annotation track (RepeatMasker output converted to BED/GFF3)
    - Protein-to-genome alignments (miniprot GFF3)