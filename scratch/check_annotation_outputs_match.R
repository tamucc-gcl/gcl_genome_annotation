library(tidyverse)
library(rtracklayer)
library(S4Vectors)

gff <- rtracklayer::import("results/<assembly_id>.annotated.gff3")
annot <- readr::read_tsv("results/<assembly_id>.functional_annotations.tsv")

mrna_ids <- gff[gff$type == "mRNA"] |>
  mcols() |>
  as.data.frame() |>
  pull(ID)

gene_ids <- gff[gff$type == "gene"] |>
  mcols() |>
  as.data.frame() |>
  pull(ID)

tibble(
  n_gff_mrna = length(mrna_ids),
  n_tsv_protein_id = n_distinct(annot$protein_id),
  n_tsv_ids_not_in_gff = sum(!annot$protein_id %in% mrna_ids),
  n_gff_ids_not_in_tsv = sum(!mrna_ids %in% annot$protein_id),
  n_gff_genes = length(gene_ids),
  n_tsv_gene_id = n_distinct(annot$gene_id),
  n_tsv_genes_not_in_gff = sum(!unique(annot$gene_id) %in% gene_ids)
)