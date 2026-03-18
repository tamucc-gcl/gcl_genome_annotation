#!/usr/bin/env Rscript
# ============================================================================
# annotation_report_plots.R
# ============================================================================
# Generate PNG plots for the genome annotation report.
#
# Usage:
#   Rscript annotation_report_plots.R \
#     --assembly_id ID \
#     --gff3 annotated.gff3 \
#     --gene_summary gene_product_summary.tsv \
#     --functional_tsv functional_annotations.tsv \
#     --repeat_tbl repeat.tbl \
#     --repeat_tbl_full buildSummary.tbl \
#     --repeat_divsum divsum \
#     --busco_txt busco.txt \
#     --trnascan_tsv trnascan.tsv \
#     --infernal_tblout infernal.tblout \
#     --telomere_gff3 telomeres.gff3 \
#     --centromere_gff3 centromere.gff3 \
#     --outdir plots/
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# ── CLI argument parsing ─────────────────────────────────────────────────────
option_list <- list(
  make_option("--assembly_id",    type = "character"),
  make_option("--gff3",           type = "character"),
  make_option("--gene_summary",   type = "character"),
  make_option("--functional_tsv", type = "character"),
  make_option("--repeat_tbl",     type = "character"),
  make_option("--repeat_tbl_full",type = "character"),
  make_option("--repeat_divsum",  type = "character"),
  make_option("--busco_txt",      type = "character"),
  make_option("--omark_summary",  type = "character"),
  make_option("--omark_detailed", type = "character"),
  make_option("--trnascan_tsv",   type = "character"),
  make_option("--infernal_tblout",type = "character"),
  make_option("--telomere_gff3",  type = "character"),
  make_option("--telomere_tsv",   type = "character"),
  make_option("--centromere_gff3",type = "character"),
  make_option("--outdir",         type = "character", default = "plots")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Null-coalescing helper ───────────────────────────────────────────────────
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ── Shared theme ─────────────────────────────────────────────────────────────
theme_report <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ── Helper: save plot ────────────────────────────────────────────────────────
save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(
    filename = file.path(opt$outdir, paste0(name, ".png")),
    plot = p, width = w, height = h, dpi = 150, bg = "white"
  )
}

# ============================================================================
# Parse GFF3
# ============================================================================
parse_gff3 <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(lines) == 0) return(data.frame())

  df <- read.delim(
    textConnection(lines), header = FALSE, sep = "\t",
    col.names = c("seqid","source","type","start","end","score","strand","phase","attributes"),
    stringsAsFactors = FALSE, quote = ""
  )
  df$length <- df$end - df$start + 1
  df
}

gff <- parse_gff3(opt$gff3)

# Split by source for convenience
braker  <- gff %>% filter(source == "BRAKER3")
trna_gff <- gff %>% filter(source == "tRNAscan-SE")
ncrna_gff <- gff %>% filter(source == "Infernal")
repeat_gff <- gff %>% filter(source == "RepeatMasker")
telo_gff  <- gff %>% filter(source == "tidk")

# ============================================================================
# SECTION A: Understanding the Annotation
# ============================================================================

# ── A1. Gene model structure distributions ───────────────────────────────────
#    Exon count per gene, gene length, CDS length, intron length
genes   <- braker %>% filter(type == "gene")
mrnas   <- braker %>% filter(type %in% c("mRNA", "transcript"))
exons   <- braker %>% filter(type == "exon")
cds     <- braker %>% filter(type == "CDS")
introns <- braker %>% filter(type == "intron")

# Extract Parent from attributes for exon counting
extract_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regmatches(attrs, regexpr(pattern, attrs, perl = TRUE))
  sub(paste0("^;?", key, "="), "", m)
}

if (nrow(exons) > 0) {
  exons$parent <- extract_attr(exons$attributes, "Parent")
  exons_per_mrna <- exons %>% count(parent, name = "n_exons")

  p_exon_dist <- ggplot(exons_per_mrna, aes(x = n_exons)) +
    geom_histogram(binwidth = 1, fill = "#3574b0", color = "white", alpha = 0.85) +
    scale_x_continuous(breaks = seq(0, max(exons_per_mrna$n_exons), by = 5)) +
    labs(
      title = "Exons per transcript",
      subtitle = paste0("Median: ", median(exons_per_mrna$n_exons),
                         ", Single-exon: ", sum(exons_per_mrna$n_exons == 1),
                         " (", round(100 * mean(exons_per_mrna$n_exons == 1), 1), "%)"),
      x = "Number of exons", y = "Count"
    ) + theme_report
  save_plot(p_exon_dist, "A1_exons_per_transcript")
}

# Gene length distribution
if (nrow(genes) > 0) {
  p_gene_len <- ggplot(genes, aes(x = length)) +
    geom_histogram(bins = 80, fill = "#3574b0", color = "white", alpha = 0.85) +
    scale_x_log10(labels = comma) +
    labs(
      title = "Gene length distribution",
      subtitle = paste0("n=", nrow(genes), ", median=", comma(median(genes$length)), " bp"),
      x = "Gene length (bp, log scale)", y = "Count"
    ) + theme_report
  save_plot(p_gene_len, "A1_gene_length_distribution")
}

# Exon length distribution
if (nrow(exons) > 0) {
  p_exon_len <- ggplot(exons, aes(x = length)) +
    geom_histogram(bins = 80, fill = "#e8913a", color = "white", alpha = 0.85) +
    scale_x_log10(labels = comma) +
    labs(
      title = "Exon length distribution",
      subtitle = paste0("n=", nrow(exons), ", median=", comma(median(exons$length)), " bp"),
      x = "Exon length (bp, log scale)", y = "Count"
    ) + theme_report
  save_plot(p_exon_len, "A1_exon_length_distribution")
}

# Intron length distribution
if (nrow(introns) > 0) {
  p_intron_len <- ggplot(introns, aes(x = length)) +
    geom_histogram(bins = 80, fill = "#6a994e", color = "white", alpha = 0.85) +
    scale_x_log10(labels = comma) +
    labs(
      title = "Intron length distribution",
      subtitle = paste0("n=", nrow(introns), ", median=", comma(median(introns$length)), " bp"),
      x = "Intron length (bp, log scale)", y = "Count"
    ) + theme_report
  save_plot(p_intron_len, "A1_intron_length_distribution")
}

# CDS length distribution
if (nrow(cds) > 0) {
  # Sum CDS length per mRNA (total coding length)
  cds$parent <- extract_attr(cds$attributes, "Parent")
  cds_per_mrna <- cds %>% group_by(parent) %>% summarise(total_cds = sum(length), .groups = "drop")

  p_cds_len <- ggplot(cds_per_mrna, aes(x = total_cds)) +
    geom_histogram(bins = 80, fill = "#a855f7", color = "white", alpha = 0.85) +
    scale_x_log10(labels = comma) +
    labs(
      title = "Total CDS length per transcript",
      subtitle = paste0("n=", nrow(cds_per_mrna), ", median=", comma(median(cds_per_mrna$total_cds)), " bp"),
      x = "CDS length (bp, log scale)", y = "Count"
    ) + theme_report
  save_plot(p_cds_len, "A1_cds_length_distribution")
}

# ── A2. Functional annotation tool coverage ──────────────────────────────────
func <- read.delim(opt$functional_tsv, stringsAsFactors = FALSE)

if (nrow(func) > 0) {
  coverage <- data.frame(
    Tool = c("InterProScan", "eggNOG", "DIAMOND", "KofamScan", "DeepSig"),
    Annotated = c(
      sum(func$ipr_signatures != "-", na.rm = TRUE),
      sum(func$eggnog_OGs != "-", na.rm = TRUE),
      sum(func$diamond_hit != "-", na.rm = TRUE),
      sum(func$kofam_KO != "-", na.rm = TRUE),
      sum(func$deepsig_type != "-" & func$deepsig_type != "Other", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  coverage$Unannotated <- nrow(func) - coverage$Annotated
  coverage$Pct <- round(100 * coverage$Annotated / nrow(func), 1)

  cov_long <- coverage %>%
    select(Tool, Annotated, Unannotated) %>%
    pivot_longer(-Tool, names_to = "Status", values_to = "Count")

  p_func_cov <- ggplot(cov_long, aes(x = reorder(Tool, -Count), y = Count, fill = Status)) +
    geom_col(position = "stack", alpha = 0.85) +
    geom_text(
      data = coverage,
      aes(x = Tool, y = Annotated, label = paste0(Pct, "%"), fill = NULL),
      vjust = -0.3, size = 3.5
    ) +
    scale_fill_manual(values = c(Annotated = "#3574b0", Unannotated = "#d1d5db")) +
    labs(
      title = "Functional annotation coverage",
      subtitle = paste0(nrow(func), " total proteins"),
      x = NULL, y = "Number of proteins"
    ) + theme_report
  save_plot(p_func_cov, "A2_functional_coverage")

  # ── Tool overlap (UpSet-style as stacked bar of annotation depth) ──────
  func$n_tools <- (
    (func$ipr_signatures != "-") +
    (func$eggnog_OGs != "-") +
    (func$diamond_hit != "-") +
    (func$kofam_KO != "-") +
    (func$deepsig_type != "-" & func$deepsig_type != "Other")
  )

  tool_depth <- func %>% count(n_tools) %>% mutate(pct = round(100 * n / sum(n), 1))

  p_tool_depth <- ggplot(tool_depth, aes(x = factor(n_tools), y = n)) +
    geom_col(fill = "#3574b0", alpha = 0.85) +
    geom_text(aes(label = paste0(pct, "%")), vjust = -0.3, size = 3.5) +
    labs(
      title = "Annotation depth per protein",
      subtitle = "Number of tools producing a hit",
      x = "Number of tools with hits", y = "Proteins"
    ) + theme_report
  save_plot(p_tool_depth, "A2_annotation_depth")
}

# ── A3. Product name composition ─────────────────────────────────────────────
gene_summary <- read.delim(opt$gene_summary, stringsAsFactors = FALSE)

if (nrow(gene_summary) > 0) {
  n_hypo <- sum(gene_summary$product == "hypothetical protein")
  n_named <- nrow(gene_summary) - n_hypo

  product_comp <- data.frame(
    Category = c("Named product", "Hypothetical protein"),
    Count = c(n_named, n_hypo)
  )

  p_product <- ggplot(product_comp, aes(x = "", y = Count, fill = Category)) +
    geom_col(width = 1, alpha = 0.85) +
    coord_polar("y") +
    scale_fill_manual(values = c("Named product" = "#3574b0", "Hypothetical protein" = "#d1d5db")) +
    labs(
      title = "Gene product assignment",
      subtitle = paste0(n_named, " named (", round(100 * n_named / nrow(gene_summary), 1),
                         "%), ", n_hypo, " hypothetical")
    ) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
  save_plot(p_product, "A3_product_composition", w = 6, h = 6)

  # Top 20 most common named products
  top_products <- gene_summary %>%
    filter(product != "hypothetical protein") %>%
    count(product, sort = TRUE) %>%
    slice_head(n = 20)

  if (nrow(top_products) > 0) {
    p_top_prod <- ggplot(top_products, aes(x = reorder(product, n), y = n)) +
      geom_col(fill = "#3574b0", alpha = 0.85) +
      coord_flip() +
      labs(
        title = "Top 20 most frequent gene products",
        x = NULL, y = "Count"
      ) + theme_report
    save_plot(p_top_prod, "A3_top_products", w = 10, h = 7)
  }
}

# ── A4. GO term distribution (top terms per ontology) ────────────────────────
if (nrow(func) > 0) {
  # Collect all GO terms from both IPR and eggNOG
  go_raw <- c(func$ipr_go, func$eggnog_GO)
  go_raw <- go_raw[go_raw != "-" & !is.na(go_raw)]
  all_go <- unlist(strsplit(go_raw, "[;,|]+"))
  all_go <- trimws(all_go)
  all_go <- all_go[grepl("^GO:", all_go)]

  if (length(all_go) > 0) {
    go_counts <- as.data.frame(table(all_go), stringsAsFactors = FALSE)
    names(go_counts) <- c("GO_term", "Count")
    go_counts <- go_counts %>% arrange(desc(Count)) %>% slice_head(n = 30)

    p_go <- ggplot(go_counts, aes(x = reorder(GO_term, Count), y = Count)) +
      geom_col(fill = "#6a994e", alpha = 0.85) +
      coord_flip() +
      labs(
        title = "Top 30 GO terms",
        subtitle = paste0(length(unique(all_go)), " unique GO terms total"),
        x = NULL, y = "Occurrences"
      ) + theme_report
    save_plot(p_go, "A4_top_go_terms", w = 9, h = 8)
  }
}

# ── A5. KEGG pathway coverage ────────────────────────────────────────────────
if (nrow(func) > 0) {
  kegg_raw <- c(func$eggnog_KEGG_path, func$kofam_KO)
  kegg_raw <- kegg_raw[kegg_raw != "-" & !is.na(kegg_raw)]
  all_kegg <- unlist(strsplit(kegg_raw, "[;,]+"))
  all_kegg <- trimws(all_kegg)
  all_kegg <- all_kegg[nchar(all_kegg) > 0]

  if (length(all_kegg) > 0) {
    kegg_counts <- as.data.frame(table(all_kegg), stringsAsFactors = FALSE)
    names(kegg_counts) <- c("KEGG_term", "Count")
    kegg_counts <- kegg_counts %>% arrange(desc(Count)) %>% slice_head(n = 30)

    p_kegg <- ggplot(kegg_counts, aes(x = reorder(KEGG_term, Count), y = Count)) +
      geom_col(fill = "#e8913a", alpha = 0.85) +
      coord_flip() +
      labs(
        title = "Top 30 KEGG terms",
        subtitle = paste0(length(unique(all_kegg)), " unique KEGG terms total"),
        x = NULL, y = "Occurrences"
      ) + theme_report
    save_plot(p_kegg, "A5_top_kegg_terms", w = 9, h = 8)
  }
}

# ── A6. InterPro domain frequency ───────────────────────────────────────────
if (nrow(func) > 0) {
  ipr_raw <- func$ipr_accessions[func$ipr_accessions != "-" & !is.na(func$ipr_accessions)]
  all_ipr <- unlist(strsplit(ipr_raw, "[;,]+"))
  all_ipr <- trimws(all_ipr)
  all_ipr <- all_ipr[grepl("^IPR", all_ipr)]

  if (length(all_ipr) > 0) {
    ipr_counts <- as.data.frame(table(all_ipr), stringsAsFactors = FALSE)
    names(ipr_counts) <- c("IPR", "Count")
    ipr_counts <- ipr_counts %>% arrange(desc(Count)) %>% slice_head(n = 30)

    p_ipr <- ggplot(ipr_counts, aes(x = reorder(IPR, Count), y = Count)) +
      geom_col(fill = "#a855f7", alpha = 0.85) +
      coord_flip() +
      labs(
        title = "Top 30 InterPro domains",
        subtitle = paste0(length(unique(all_ipr)), " unique IPR accessions total"),
        x = NULL, y = "Occurrences"
      ) + theme_report
    save_plot(p_ipr, "A6_top_interpro_domains", w = 9, h = 8)
  }
}

# ── A7. tRNA isotype distribution ────────────────────────────────────────────
if (nrow(trna_gff) > 0) {
  trna_features <- trna_gff %>% filter(type == "tRNA")
  trna_features$isotype <- sub(".*isotype=([^;]+).*", "\\1", trna_features$attributes)
  trna_features$isotype[!grepl("isotype=", trna_features$attributes)] <- "Unknown"

  iso_counts <- trna_features %>% count(isotype, sort = TRUE)

  p_trna <- ggplot(iso_counts, aes(x = reorder(isotype, n), y = n)) +
    geom_col(fill = "#e05c5c", alpha = 0.85) +
    coord_flip() +
    labs(
      title = "tRNA isotype distribution",
      subtitle = paste0(nrow(trna_features), " tRNAs total"),
      x = "Isotype", y = "Count"
    ) + theme_report
  save_plot(p_trna, "A7_trna_isotypes", w = 7, h = 6)
}

# ── A8. ncRNA Rfam family distribution ───────────────────────────────────────
if (nrow(ncrna_gff) > 0) {
  ncrna_gff$rfam_name <- sub(".*Name=([^;]+).*", "\\1", ncrna_gff$attributes)
  rfam_counts <- ncrna_gff %>% count(rfam_name, sort = TRUE) %>% slice_head(n = 25)

  p_ncrna <- ggplot(rfam_counts, aes(x = reorder(rfam_name, n), y = n)) +
    geom_col(fill = "#f59e42", alpha = 0.85) +
    coord_flip() +
    labs(
      title = "Top 25 Rfam ncRNA families",
      subtitle = paste0(nrow(ncrna_gff), " ncRNA features total"),
      x = "Rfam family", y = "Count"
    ) + theme_report
  save_plot(p_ncrna, "A8_ncrna_rfam_families", w = 8, h = 7)
}

# ── A9. Telomere and centromere positional summary ───────────────────────────
if (nrow(telo_gff) > 0) {
  telo_gff$position <- sub(".*position=([^;]+).*", "\\1", telo_gff$attributes)
  telo_gff$position[!grepl("position=", telo_gff$attributes)] <- "unknown"

  telo_pos <- telo_gff %>% count(position)

  p_telo_pos <- ggplot(telo_pos, aes(x = reorder(position, -n), y = n, fill = position)) +
    geom_col(alpha = 0.85) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Telomere region classification",
      subtitle = paste0(nrow(telo_gff), " telomeric regions"),
      x = "Position class", y = "Count"
    ) + theme_report + theme(legend.position = "none")
  save_plot(p_telo_pos, "A9_telomere_positions", w = 6, h = 5)
}

# ── A9b. Telomere density along scaffolds (from tidk TSV) ───────────────────
#   tidk search output columns: id, window, forward_repeat_number,
#   reverse_repeat_number, telomeric_repeat
telo_tsv <- tryCatch(
  read.delim(opt$telomere_tsv, stringsAsFactors = FALSE),
  error = function(e) data.frame()
)

if (nrow(telo_tsv) > 0) {
  # Standardize column names (tidk uses varying names across versions)
  names(telo_tsv) <- gsub("\\.", "_", tolower(names(telo_tsv)))
  if (!"id" %in% names(telo_tsv) && length(names(telo_tsv)) >= 1) {
    names(telo_tsv)[1] <- "id"
  }

  # Identify repeat count columns
  fwd_col <- grep("forward", names(telo_tsv), value = TRUE)[1]
  rev_col <- grep("reverse", names(telo_tsv), value = TRUE)[1]
  win_col <- grep("window", names(telo_tsv), value = TRUE)[1]

  if (!is.na(fwd_col) && !is.na(rev_col) && !is.na(win_col)) {
    telo_tsv$total_repeats <- telo_tsv[[fwd_col]] + telo_tsv[[rev_col]]
    telo_tsv$position_mb <- telo_tsv[[win_col]] / 1e6

    # Compute scaffold lengths and select top N scaffolds by length
    scaff_max <- telo_tsv %>%
      group_by(id) %>%
      summarise(max_pos = max(.data[[win_col]]), .groups = "drop") %>%
      arrange(desc(max_pos))

    # Show top 24 scaffolds (fits a 6x4 facet grid)
    top_scaffolds <- scaff_max %>% slice_head(n = 24) %>% pull(id)
    telo_plot_data <- telo_tsv %>% filter(id %in% top_scaffolds)

    # Order scaffolds by length
    telo_plot_data$id <- factor(telo_plot_data$id, levels = top_scaffolds)

    if (nrow(telo_plot_data) > 0) {
      p_telo_density <- ggplot(telo_plot_data, aes(x = position_mb, y = total_repeats)) +
        geom_area(fill = "#e05c5c", alpha = 0.5) +
        geom_line(color = "#e05c5c", linewidth = 0.3) +
        facet_wrap(~ id, scales = "free_x", ncol = 4) +
        labs(
          title = "Telomere repeat density along scaffolds",
          subtitle = paste0("Top ", length(top_scaffolds),
                            " scaffolds by length (tidk windowed counts)"),
          x = "Position (Mb)", y = "Repeats per window"
        ) +
        theme_report +
        theme(
          strip.text = element_text(size = 7),
          axis.text = element_text(size = 6),
          panel.spacing = unit(0.3, "lines")
        )

      n_rows <- ceiling(length(top_scaffolds) / 4)
      save_plot(p_telo_density, "A9b_telomere_density",
                w = 12, h = max(4, n_rows * 2.2))
    }
  }
}

# ── A10. Per-scaffold annotation density ─────────────────────────────────────
if (nrow(genes) > 0) {
  scaff_genes <- genes %>% count(seqid, name = "n_genes") %>% arrange(desc(n_genes))

  # Top 30 scaffolds by gene count
  top_scaff <- scaff_genes %>% slice_head(n = 30)

  p_scaff_density <- ggplot(top_scaff, aes(x = reorder(seqid, n_genes), y = n_genes)) +
    geom_col(fill = "#3574b0", alpha = 0.85) +
    coord_flip() +
    labs(
      title = "Genes per scaffold (top 30)",
      subtitle = paste0(n_distinct(genes$seqid), " scaffolds with genes"),
      x = NULL, y = "Gene count"
    ) + theme_report
  save_plot(p_scaff_density, "A10_genes_per_scaffold", w = 9, h = 7)
}

# ============================================================================
# SECTION B: Quality Assessment
# ============================================================================

# ── B1. BUSCO completeness bar ───────────────────────────────────────────────
busco_lines <- readLines(opt$busco_txt)
busco_summary_line <- grep("^\\s*C:", busco_lines, value = TRUE)

if (length(busco_summary_line) > 0) {
  # Parse C:XX.X%[S:XX.X%,D:XX.X%],F:XX.X%,M:XX.X%,n:NNNN
  bs <- busco_summary_line[1]
  complete_s <- as.numeric(sub(".*S:([0-9.]+)%.*", "\\1", bs))
  complete_d <- as.numeric(sub(".*D:([0-9.]+)%.*", "\\1", bs))
  fragmented <- as.numeric(sub(".*F:([0-9.]+)%.*", "\\1", bs))
  missing    <- as.numeric(sub(".*M:([0-9.]+)%.*", "\\1", bs))

  busco_df <- data.frame(
    Category = factor(
      c("Complete (single)", "Complete (duplicated)", "Fragmented", "Missing"),
      levels = c("Missing", "Fragmented", "Complete (duplicated)", "Complete (single)")
    ),
    Pct = c(complete_s, complete_d, fragmented, missing)
  )

  p_busco <- ggplot(busco_df, aes(x = "BUSCO", y = Pct, fill = Category)) +
    geom_col(width = 0.6, alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Complete (single)" = "#3574b0",
      "Complete (duplicated)" = "#6ab0e3",
      "Fragmented" = "#e8913a",
      "Missing" = "#d1d5db"
    )) +
    labs(
      title = "BUSCO completeness",
      subtitle = paste0("C:", complete_s + complete_d, "% [S:", complete_s,
                         "%, D:", complete_d, "%], F:", fragmented, "%, M:", missing, "%"),
      x = NULL, y = "Percentage"
    ) +
    theme_report + theme(axis.text.y = element_blank())
  save_plot(p_busco, "B1_busco_completeness", w = 8, h = 3.5)
}

# ── B1b. OMArk completeness and consistency ──────────────────────────────────
#   OMArk summary format (compact lines):
#     S:8681,D:3569[U:3555,E:14],M:2383           (counts)
#     S:59.32%,D:24.39%[U:24.29%,E:0.10%],M:16.29% (percentages)
#     A:27432[P:5164,F:3715],I:3001[...],C:0[...],U:13696  (counts)
#     A:62.16%[P:11.70%,...],I:6.80%[...],C:0.00%[...],U:31.04% (percentages)
omark_lines <- tryCatch(readLines(opt$omark_summary), error = function(e) character(0))

if (length(omark_lines) > 0) {

  # Helper: extract number after a key like "S:" from a compact line
  extract_omark_pct <- function(lines, key) {
    for (line in lines) {
      pat <- paste0(key, "([0-9]+\\.?[0-9]*)%")
      m <- regmatches(line, regexpr(pat, line, perl = TRUE))
      if (length(m) > 0) {
        val <- sub(paste0("^", key), "", m[1])
        val <- sub("%$", "", val)
        return(as.numeric(val))
      }
    }
    return(0)
  }

  # ── Completeness: S / D(U+E) / M ──────────────────────────────────────
  s_pct  <- extract_omark_pct(omark_lines, "S:")
  d_pct  <- extract_omark_pct(omark_lines, "D:")
  du_pct <- extract_omark_pct(omark_lines, "U:")
  de_pct <- extract_omark_pct(omark_lines, "E:")
  m_pct  <- extract_omark_pct(omark_lines, "M:")

  if (s_pct + d_pct + m_pct > 0) {
    comp_df <- data.frame(
      Category = factor(
        c("Single", "Duplicated (unexpected)", "Duplicated (expected)", "Missing"),
        levels = rev(c("Single", "Duplicated (unexpected)", "Duplicated (expected)", "Missing"))
      ),
      Pct = c(s_pct, du_pct, de_pct, m_pct)
    ) %>% filter(Pct > 0)

    comp_subtitle <- paste0(
      "Clade: ",
      sub(".*selected clade was\\s+", "",
          grep("selected clade", omark_lines, value = TRUE, ignore.case = TRUE)[1] %||% ""),
      "  |  S:", s_pct, "%, D:", d_pct, "% [U:", du_pct, "%, E:", de_pct, "%], M:", m_pct, "%"
    )

    p_omark_comp <- ggplot(comp_df, aes(x = "Completeness", y = Pct, fill = Category)) +
      geom_col(width = 0.6, alpha = 0.85) +
      geom_text(aes(label = ifelse(Pct >= 2, paste0(Pct, "%"), "")),
                position = position_stack(vjust = 0.5), size = 3.5) +
      coord_flip() +
      scale_fill_manual(values = c(
        "Single" = "#3574b0",
        "Duplicated (unexpected)" = "#e8913a",
        "Duplicated (expected)" = "#6ab0e3",
        "Missing" = "#d1d5db"
      )) +
      labs(title = "OMArk — Conserved gene completeness",
           subtitle = comp_subtitle,
           x = NULL, y = "Percentage") +
      theme_report + theme(axis.text.y = element_blank())
    save_plot(p_omark_comp, "B1b_omark_completeness", w = 9, h = 3.5)
  }

  # ── Consistency: A / I / C / U ─────────────────────────────────────────
  a_pct  <- extract_omark_pct(omark_lines, "A:")
  i_pct  <- extract_omark_pct(omark_lines, "I:")
  c_pct  <- extract_omark_pct(omark_lines, "C:")
  u_pct  <- extract_omark_pct(omark_lines, "U:")

  if (a_pct + i_pct + c_pct + u_pct > 0) {
    # Also extract sub-categories for the consistent portion
    ap_pct <- extract_omark_pct(omark_lines, "P:")
    af_pct <- extract_omark_pct(omark_lines, "F:")

    cons_df <- data.frame(
      Category = factor(
        c("Consistent", "Inconsistent", "Contaminant", "Unknown"),
        levels = rev(c("Consistent", "Inconsistent", "Contaminant", "Unknown"))
      ),
      Pct = c(a_pct, i_pct, c_pct, u_pct)
    ) %>% filter(Pct > 0)

    # Extract total protein count from summary
    n_prot_line <- grep("there are.*proteins", omark_lines, value = TRUE, ignore.case = TRUE)[1] %||% ""
    n_prot_str <- sub(".*there are\\s+", "", n_prot_line)
    n_prot_str <- sub("\\s+proteins.*", "", n_prot_str)

    cons_subtitle <- paste0(
      n_prot_str, " proteins  |  A:", a_pct, "%, I:", i_pct,
      "%, C:", c_pct, "%, U:", u_pct, "%"
    )

    p_omark_cons <- ggplot(cons_df, aes(x = "Consistency", y = Pct, fill = Category)) +
      geom_col(width = 0.6, alpha = 0.85) +
      geom_text(aes(label = ifelse(Pct >= 2, paste0(Pct, "%"), "")),
                position = position_stack(vjust = 0.5), size = 3.5) +
      coord_flip() +
      scale_fill_manual(values = c(
        "Consistent" = "#6a994e",
        "Inconsistent" = "#e8913a",
        "Contaminant" = "#dc2626",
        "Unknown" = "#9ca3af"
      )) +
      labs(title = "OMArk — Proteome consistency",
           subtitle = cons_subtitle,
           x = NULL, y = "Percentage") +
      theme_report + theme(axis.text.y = element_blank())
    save_plot(p_omark_cons, "B1c_omark_consistency", w = 9, h = 3.5)

    # ── Consistency sub-breakdown (partial/fragmented within consistent) ──
    if (a_pct > 0) {
      a_full_pct <- a_pct - ap_pct - af_pct
      # Same for inconsistent if present
      # P: and F: appear multiple times — need context-aware parsing
      # For now use the detailed file for the sub-breakdown

      sub_df <- data.frame(
        Category = factor(
          c("Full-length", "Partial hit", "Fragmented"),
          levels = c("Full-length", "Partial hit", "Fragmented")
        ),
        Pct = c(max(0, a_full_pct), ap_pct, af_pct)
      ) %>% filter(Pct > 0)

      if (nrow(sub_df) > 1) {
        p_omark_sub <- ggplot(sub_df, aes(x = Category, y = Pct, fill = Category)) +
          geom_col(width = 0.6, alpha = 0.85) +
          geom_text(aes(label = paste0(Pct, "%")), vjust = -0.3, size = 3.5) +
          scale_fill_manual(values = c(
            "Full-length" = "#6a994e",
            "Partial hit" = "#eab308",
            "Fragmented" = "#dc2626"
          )) +
          labs(
            title = "OMArk — Consistent gene quality",
            subtitle = paste0("Breakdown of ", a_pct, "% consistent proteins"),
            x = NULL, y = "% of proteome"
          ) + theme_report + theme(legend.position = "none")
        save_plot(p_omark_sub, "B1d_omark_consistent_quality", w = 6, h = 5)
      }
    }
  }
}

# ── B2. Short intron / short exon QC flags ───────────────────────────────────
if (nrow(introns) > 0) {
  intron_qc <- data.frame(
    Category = c("< 10 bp (suspicious)", "10-30 bp (short)", "31-100 bp", "> 100 bp"),
    Count = c(
      sum(introns$length < 10),
      sum(introns$length >= 10 & introns$length <= 30),
      sum(introns$length > 30 & introns$length <= 100),
      sum(introns$length > 100)
    )
  )
  intron_qc$Category <- factor(intron_qc$Category, levels = intron_qc$Category)

  p_intron_qc <- ggplot(intron_qc, aes(x = Category, y = Count, fill = Category)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c(
      "< 10 bp (suspicious)" = "#dc2626",
      "10-30 bp (short)" = "#e8913a",
      "31-100 bp" = "#eab308",
      "> 100 bp" = "#6a994e"
    )) +
    labs(
      title = "Intron length QC",
      subtitle = "Introns < 10 bp are likely artifacts; canonical minimum ~30 bp",
      x = NULL, y = "Count"
    ) + theme_report + theme(legend.position = "none")
  save_plot(p_intron_qc, "B2_intron_length_qc", w = 7, h = 5)
}

# ── B3. Gene length vs exon count (diagnostic scatter) ──────────────────────
if (nrow(genes) > 0 && exists("exons_per_mrna")) {
  genes$gene_id <- extract_attr(genes$attributes, "ID")

  mrnas_info <- mrnas %>%
    mutate(mrna_id = extract_attr(attributes, "ID"),
           parent_gene = extract_attr(attributes, "Parent"))

  genes_lookup <- genes %>% select(gene_id, gene_length = length)

  gene_vs_exon <- mrnas_info %>%
    left_join(exons_per_mrna, by = c("mrna_id" = "parent")) %>%
    left_join(genes_lookup, by = c("parent_gene" = "gene_id"))

  if (nrow(gene_vs_exon) > 0 && !all(is.na(gene_vs_exon$n_exons))) {
    p_scatter <- ggplot(gene_vs_exon, aes(x = n_exons, y = gene_length)) +
      geom_point(alpha = 0.15, size = 0.8, color = "#3574b0") +
      scale_y_log10(labels = comma) +
      labs(
        title = "Gene length vs exon count",
        subtitle = "Outliers may indicate gene model errors",
        x = "Exons per transcript", y = "Gene length (bp, log scale)"
      ) + theme_report
    save_plot(p_scatter, "B3_gene_length_vs_exons", w = 7, h = 6)
  }
}

# ── B4. Functional annotation agreement matrix ──────────────────────────────
#   How often do tools agree on annotating the same protein?
if (nrow(func) > 0) {
  tools <- c("ipr_signatures", "eggnog_OGs", "diamond_hit", "kofam_KO")
  tool_names <- c("InterProScan", "eggNOG", "DIAMOND", "KofamScan")
  n_tools <- length(tools)

  overlap_mat <- matrix(0, nrow = n_tools, ncol = n_tools,
                        dimnames = list(tool_names, tool_names))

  for (i in seq_len(n_tools)) {
    for (j in seq_len(n_tools)) {
      hit_i <- func[[tools[i]]] != "-"
      hit_j <- func[[tools[j]]] != "-"
      overlap_mat[i, j] <- sum(hit_i & hit_j)
    }
  }

  overlap_df <- as.data.frame(as.table(overlap_mat))
  names(overlap_df) <- c("Tool1", "Tool2", "Count")

  p_overlap <- ggplot(overlap_df, aes(x = Tool1, y = Tool2, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = comma(Count)), size = 3.5) +
    scale_fill_gradient(low = "white", high = "#3574b0") +
    labs(
      title = "Functional tool overlap",
      subtitle = "Pairwise count of proteins annotated by both tools",
      x = NULL, y = NULL
    ) + theme_report +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  save_plot(p_overlap, "B4_tool_overlap_matrix", w = 6, h = 6)
}

# ── B5. Repeat composition summary ──────────────────────────────────────────
if (nrow(repeat_gff) > 0) {
  # Extract repeat class from Target attribute or type field
  repeat_gff$repeat_class <- sub(".*Target=([^ ]+).*", "\\1", repeat_gff$attributes)
  repeat_gff$repeat_class <- sub("#.*", "", repeat_gff$repeat_class)

  # Simplify classes
  repeat_gff$main_class <- sub("/.*", "", repeat_gff$repeat_class)

  repeat_summary <- repeat_gff %>%
    group_by(main_class) %>%
    summarise(total_bp = sum(length), count = n(), .groups = "drop") %>%
    arrange(desc(total_bp))

  p_repeat <- ggplot(repeat_summary %>% slice_head(n = 12),
                     aes(x = reorder(main_class, total_bp), y = total_bp / 1e6)) +
    geom_col(fill = "#8b5e3c", alpha = 0.85) +
    coord_flip() +
    labs(
      title = "Repeat composition by class",
      subtitle = paste0(comma(sum(repeat_summary$total_bp)), " bp total repeat content"),
      x = NULL, y = "Total length (Mbp)"
    ) + theme_report
  save_plot(p_repeat, "B5_repeat_composition", w = 8, h = 6)
}

# ============================================================================
# Summary message
# ============================================================================
plot_files <- list.files(opt$outdir, pattern = "\\.png$")
message(sprintf("Generated %d plots in %s/", length(plot_files), opt$outdir))