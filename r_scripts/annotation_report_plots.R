#!/usr/bin/env Rscript
# ============================================================================
# annotation_report_plots.R
# ============================================================================
# Generate PNG plots for a unified multi-assembly annotation report.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(treemapify)
  library(GO.db)
  library(AnnotationDbi)
})

# Restore dplyr::select after AnnotationDbi masks it
select <- dplyr::select

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--metrics",          type = "character"),
  make_option("--gff3s",            type = "character"),
  make_option("--gene_summaries",   type = "character"),
  make_option("--functional_tsvs",  type = "character"),
  make_option("--busco_txts",       type = "character"),
  make_option("--omark_summaries",  type = "character"),
  make_option("--trnascan_tsvs",    type = "character"),
  make_option("--infernal_tblouts", type = "character"),
  make_option("--telomere_gff3s",   type = "character"),
  make_option("--telomere_tsvs",    type = "character"),
  make_option("--centromere_gff3s", type = "character"),
  make_option("--outdir",           type = "character", default = "plots")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# ── Helpers ──────────────────────────────────────────────────────────────────
theme_report <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

theme_compare <- theme_report +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))

save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(opt$outdir, paste0(name, ".png")),
         plot = p, width = w, height = h, dpi = 150, bg = "white")
}

should_log <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) return(FALSE)
  max(x) / min(x) > 100
}

smart_y <- function(p, values) {
  if (should_log(values)) p + scale_y_log10(labels = comma)
  else p + scale_y_continuous(labels = comma)
}

smart_x <- function(p, values) {
  if (should_log(values)) p + scale_x_log10(labels = comma)
  else p + scale_x_continuous(labels = comma)
}

smart_y_bar <- function(p, values) {
  if (should_log(values)) p + scale_y_sqrt(labels = comma)
  else p + scale_y_continuous(labels = comma)
}

split_files <- function(arg) {
  if (is.null(arg) || arg == "") return(character(0))
  strsplit(arg, "\\s+")[[1]]
}

infer_id <- function(path, suffixes) {
  base <- basename(path)
  for (s in suffixes) {
    if (endsWith(base, s)) return(sub(paste0(gsub("\\.", "\\\\.", s), "$"), "", base))
  }
  tools::file_path_sans_ext(base)
}

parse_gff3 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(lines) == 0) return(data.frame())
  df <- read.delim(textConnection(lines), header = FALSE, sep = "\t",
    col.names = c("seqid","source","type","start","end","score","strand","phase","attributes"),
    stringsAsFactors = FALSE, quote = "")
  df$length <- df$end - df$start + 1
  df
}

extract_attr <- function(attrs, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regmatches(attrs, regexpr(pattern, attrs, perl = TRUE))
  sub(paste0("^;?", key, "="), "", m)
}

# ── Read metrics TSV ─────────────────────────────────────────────────────────
metrics <- read.delim(opt$metrics, stringsAsFactors = FALSE)
if (nrow(metrics) == 0) { message("No assemblies found."); quit(save = "no") }

metrics$label <- ifelse(nchar(metrics$assembly_id) > 25,
                        paste0(substr(metrics$assembly_id, 1, 22), "..."),
                        metrics$assembly_id)

assembly_ids <- metrics$assembly_id

# ── Split file lists ─────────────────────────────────────────────────────────
gff3_files         <- split_files(opt$gff3s)
gene_summary_files <- split_files(opt$gene_summaries)
functional_files   <- split_files(opt$functional_tsvs)
telomere_tsv_files <- split_files(opt$telomere_tsvs)
telomere_gff3_files <- split_files(opt$telomere_gff3s)

# ============================================================================
# PER-ASSEMBLY PLOTS
# ============================================================================

for (gff_path in gff3_files) {
  aid <- infer_id(gff_path, c(".annotated.gff3"))
  gff <- tryCatch(parse_gff3(gff_path), error = function(e) data.frame())
  if (nrow(gff) == 0) next

  braker  <- gff %>% filter(source == "BRAKER3")
  genes   <- braker %>% filter(type == "gene")
  exons   <- braker %>% filter(type == "exon")
  introns <- braker %>% filter(type == "intron")

  # ── Gene length distribution ───────────────────────────────────────
  if (nrow(genes) > 0) {
    p <- ggplot(genes, aes(x = length)) +
      geom_histogram(bins = 80, fill = "#3574b0", color = "white", alpha = 0.85) +
      labs(title = paste0(aid, " — Gene length distribution"),
           subtitle = paste0("n=", comma(nrow(genes)),
                             ", median=", comma(median(genes$length)), " bp"),
           x = "Gene length (bp)", y = "Count") + theme_report
    p <- smart_x(p, genes$length)
    p <- p + scale_y_continuous(labels = comma)
    save_plot(p, paste0(aid, "_gene_length_dist"))
  }

  # ── Exons per transcript ───────────────────────────────────────────
  if (nrow(exons) > 0) {
    exons$parent <- extract_attr(exons$attributes, "Parent")
    epm <- exons %>% count(parent, name = "n_exons")
    cap <- 50
    n_over <- sum(epm$n_exons > cap)
    over_note <- if (n_over > 0) paste0(", ", comma(n_over), " >", cap, " not shown") else ""

    p <- ggplot(epm %>% filter(n_exons <= cap), aes(x = n_exons)) +
      geom_histogram(binwidth = 1, fill = "#3574b0", color = "white", alpha = 0.85) +
      scale_x_continuous(breaks = seq(0, cap, by = 5)) +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Exons per transcript"),
           subtitle = paste0("Median: ", median(epm$n_exons),
                             ", single-exon: ", comma(sum(epm$n_exons == 1)),
                             " (", round(100 * mean(epm$n_exons == 1), 1), "%)",
                             over_note),
           x = "Exon count", y = "Transcripts") + theme_report
    save_plot(p, paste0(aid, "_exons_per_tx"))
  }

  # ── Intron length QC ───────────────────────────────────────────────
  if (nrow(introns) > 0) {
    iq <- data.frame(
      Category = factor(c("< 10 bp", "10–30 bp", "31–100 bp", "> 100 bp"),
                        levels = c("< 10 bp", "10–30 bp", "31–100 bp", "> 100 bp")),
      Count = c(sum(introns$length < 10),
                sum(introns$length >= 10 & introns$length <= 30),
                sum(introns$length > 30 & introns$length <= 100),
                sum(introns$length > 100)))
    p <- ggplot(iq, aes(x = Category, y = Count, fill = Category)) +
      geom_col(alpha = 0.85) +
      scale_fill_manual(values = c("#dc2626","#e8913a","#eab308","#6a994e")) +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Intron length QC"),
           subtitle = "Introns < 10 bp likely artifacts",
           x = NULL, y = "Count") +
      theme_report + theme(legend.position = "none")
    save_plot(p, paste0(aid, "_intron_qc"))
  }

  # ── tRNA isotypes ──────────────────────────────────────────────────
  trna_gff <- gff %>% filter(source == "tRNAscan-SE", type == "tRNA")
  if (nrow(trna_gff) > 0) {
    trna_gff$isotype <- sub(".*isotype=([^;]+).*", "\\1", trna_gff$attributes)
    trna_gff$isotype[!grepl("isotype=", trna_gff$attributes)] <- "Unknown"
    iso <- trna_gff %>% count(isotype, sort = TRUE)
    p <- ggplot(iso, aes(x = reorder(isotype, n), y = n)) +
      geom_col(fill = "#e05c5c", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — tRNA isotypes"),
           subtitle = paste0(nrow(trna_gff), " tRNAs"),
           x = NULL, y = "Count") + theme_report
    save_plot(p, paste0(aid, "_trna_isotypes"), w = 7, h = 6)
  }

  # ── ncRNA Rfam families ────────────────────────────────────────────
  ncrna <- gff %>% filter(source == "Infernal")
  if (nrow(ncrna) > 0) {
    ncrna$rfam <- sub(".*Name=([^;]+).*", "\\1", ncrna$attributes)
    rfam_ct <- ncrna %>% count(rfam, sort = TRUE) %>% slice_head(n = 25)
    p <- ggplot(rfam_ct, aes(x = reorder(rfam, n), y = n)) +
      geom_col(fill = "#f59e42", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Top Rfam families"),
           subtitle = paste0(nrow(ncrna), " ncRNA features"),
           x = NULL, y = "Count") + theme_report
    save_plot(p, paste0(aid, "_ncrna_rfam"), w = 8, h = 7)
  }

  # ── Genes per scaffold ────────────────────────────────────────────
  if (nrow(genes) > 0) {
    sg <- genes %>% count(seqid, name = "n_genes") %>%
      arrange(desc(n_genes)) %>% slice_head(n = 30)
    p <- ggplot(sg, aes(x = reorder(seqid, n_genes), y = n_genes)) +
      geom_col(fill = "#3574b0", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Genes per scaffold (top 30)"),
           x = NULL, y = "Gene count") + theme_report
    save_plot(p, paste0(aid, "_genes_per_scaffold"), w = 9, h = 7)
  }

  # ── Per-scaffold repeat masking % ──────────────────────────────────
  repeats <- gff %>% filter(source == "RepeatMasker")
  if (nrow(repeats) > 0) {
    masked_bp <- repeats %>%
      group_by(seqid) %>%
      summarise(masked_bp = sum(length), .groups = "drop")

    scaff_lens <- gff %>%
      group_by(seqid) %>%
      summarise(scaff_len = max(end), .groups = "drop")

    sm <- masked_bp %>%
      inner_join(scaff_lens, by = "seqid") %>%
      mutate(pct_masked = pmin(100, 100 * masked_bp / scaff_len)) %>%
      filter(scaff_len >= 100000) %>%
      arrange(desc(scaff_len))

    if (nrow(sm) > 0) {
      sm$label <- ifelse(nchar(sm$seqid) > 25,
                         paste0(substr(sm$seqid, 1, 22), "..."), sm$seqid)
      sm$label <- factor(sm$label, levels = rev(sm$label))
      plot_h <- max(5, min(20, nrow(sm) * 0.3 + 2))

      p <- ggplot(sm, aes(x = label, y = pct_masked)) +
        geom_col(fill = "#8b5e3c", alpha = 0.85) +
        geom_text(aes(label = paste0(round(pct_masked, 1), "%")),
                  hjust = -0.1, size = 2.5) +
        coord_flip() +
        scale_y_continuous(limits = c(0, 105),
                           labels = function(x) paste0(x, "%")) +
        labs(title = paste0(aid, " — Repeat masking per scaffold"),
             subtitle = paste0(nrow(sm), " scaffolds >= 100 kb"),
             x = NULL, y = "% scaffold masked") + theme_report
      save_plot(p, paste0(aid, "_repeat_masking_per_scaffold"),
                w = 10, h = plot_h)
    }
  }
}

# ── Per-assembly: functional annotation ──────────────────────────────────────
for (func_path in functional_files) {
  aid <- infer_id(func_path, c(".functional_annotations.tsv"))
  func <- tryCatch(read.delim(func_path, stringsAsFactors = FALSE),
                   error = function(e) data.frame())
  if (nrow(func) == 0) next

  # ── AED distribution ───────────────────────────────────────────────
  if ("aed" %in% names(func)) {
    aed_df <- func %>%
      filter(!is.na(aed), aed != ".", aed != "-") %>%
      mutate(aed = as.numeric(aed)) %>%
      filter(!is.na(aed))

    if (nrow(aed_df) > 0) {
      # Cumulative distribution (the standard MAKER-style figure)
      p <- ggplot(aed_df, aes(x = aed)) +
        stat_ecdf(geom = "step", linewidth = 1, colour = "#2563eb") +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(labels = percent) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
        annotate("text", x = 0.52, y = 0.1, label = "AED = 0.5",
                 color = "grey50", hjust = 0, size = 3) +
        labs(title = paste0(aid, " — Cumulative AED Distribution"),
             subtitle = paste0("n=", comma(nrow(aed_df)),
                               ", median=", round(median(aed_df$aed), 3),
                               ", <0.5: ", round(100 * mean(aed_df$aed < 0.5), 1), "%"),
             x = "Annotation Edit Distance (AED)",
             y = "Cumulative fraction of transcripts") + theme_report
      save_plot(p, paste0(aid, "_aed_cumulative"))

      # Histogram
      p2 <- ggplot(aed_df, aes(x = aed)) +
        geom_histogram(binwidth = 0.05, fill = "#2563eb",
                       colour = "white", alpha = 0.85) +
        scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(labels = comma) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
        labs(title = paste0(aid, " — AED Distribution"),
             subtitle = paste0("n=", comma(nrow(aed_df)),
                               ", AED=1 (no evidence): ",
                               comma(sum(aed_df$aed >= 0.999))),
             x = "AED", y = "Number of transcripts") + theme_report
      save_plot(p2, paste0(aid, "_aed_histogram"))
    }
  }

  # Tool coverage bar
  cov <- data.frame(
    Tool = c("InterProScan", "eggNOG", "DIAMOND", "KofamScan", "DeepSig"),
    Annotated = c(
      sum(func$ipr_signatures != "-", na.rm = TRUE),
      sum(func$eggnog_OGs != "-", na.rm = TRUE),
      sum(func$diamond_hit != "-", na.rm = TRUE),
      sum(func$kofam_KO != "-", na.rm = TRUE),
      sum(func$deepsig_type != "-" & func$deepsig_type != "Other", na.rm = TRUE)
    ), stringsAsFactors = FALSE)
  cov$Pct <- round(100 * cov$Annotated / nrow(func), 1)

  p <- ggplot(cov, aes(x = reorder(Tool, -Annotated), y = Annotated)) +
    geom_col(fill = "#3574b0", alpha = 0.85) +
    geom_text(aes(label = paste0(Pct, "%")), vjust = -0.3, size = 3.5) +
    scale_y_continuous(labels = comma) +
    labs(title = paste0(aid, " — Functional tool coverage"),
         subtitle = paste0(comma(nrow(func)), " proteins"),
         x = NULL, y = "Proteins annotated") + theme_report
  save_plot(p, paste0(aid, "_functional_coverage"))

  # ── Top GO terms ───────────────────────────────────────────────────
  go_raw <- c(func$ipr_go, func$eggnog_GO)
  go_raw <- go_raw[go_raw != "-" & !is.na(go_raw)]
  all_go <- unlist(strsplit(go_raw, "[;,|]+"))
  all_go <- trimws(sub("\\(.*\\)$", "", trimws(all_go)))
  all_go <- all_go[grepl("^GO:\\d+$", all_go)]

  if (length(all_go) > 0) {
    go_ct <- as.data.frame(table(all_go), stringsAsFactors = FALSE)
    names(go_ct) <- c("GOID", "Count")
    go_ct <- go_ct %>% arrange(desc(Count)) %>% slice_head(n = 30)

    go_names <- tryCatch(
      AnnotationDbi::select(GO.db, keys = go_ct$GOID,
                            columns = "TERM", keytype = "GOID"),
      error = function(e) data.frame(GOID = character(), TERM = character()))

    go_ct <- go_ct %>%
      left_join(go_names, by = "GOID") %>%
      mutate(label = ifelse(
        !is.na(TERM) & TERM != "",
        paste0(GOID, ": ", ifelse(nchar(TERM) > 45,
               paste0(substr(TERM, 1, 42), "..."), TERM)),
        GOID))

    p <- ggplot(go_ct, aes(x = reorder(label, Count), y = Count)) +
      geom_col(fill = "#6a994e", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Top GO terms"),
           subtitle = paste0(length(unique(all_go)), " unique"),
           x = NULL, y = "Occurrences") + theme_report
    save_plot(p, paste0(aid, "_top_go"), w = 11, h = 9)

    # GO treemap
    go_all_ct <- as.data.frame(table(all_go), stringsAsFactors = FALSE)
    names(go_all_ct) <- c("GOID", "Count")
    go_annot <- tryCatch(
      AnnotationDbi::select(GO.db, keys = go_all_ct$GOID,
                            columns = c("TERM","ONTOLOGY"), keytype = "GOID"),
      error = function(e) data.frame())

    if (nrow(go_annot) > 0) {
      go_tree <- go_all_ct %>%
        inner_join(go_annot, by = "GOID") %>%
        filter(!is.na(ONTOLOGY), !is.na(TERM)) %>%
        mutate(Ontology = recode(ONTOLOGY,
          BP = "Biological Process",
          MF = "Molecular Function",
          CC = "Cellular Component"))

      go_top <- go_tree %>%
        group_by(Ontology) %>%
        slice_max(order_by = Count, n = 15, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(label = ifelse(nchar(TERM) > 40,
               paste0(substr(TERM, 1, 37), "..."), TERM))

      if (nrow(go_top) > 0) {
        p <- ggplot(go_top,
                    aes(area = Count, fill = Ontology,
                        label = label, subgroup = Ontology)) +
          geom_treemap(alpha = 0.85) +
          geom_treemap_subgroup_border(color = "white", size = 2) +
          geom_treemap_text(colour = "white", place = "centre",
                            grow = FALSE, reflow = TRUE, size = 9) +
          scale_fill_manual(values = c(
            "Biological Process" = "#3574b0",
            "Molecular Function" = "#e8913a",
            "Cellular Component" = "#6a994e")) +
          labs(title = paste0(aid, " — GO treemap"),
               subtitle = paste0("Top 15/ontology, ",
                                 length(unique(all_go)), " total")) +
          theme_report
        save_plot(p, paste0(aid, "_go_treemap"), w = 12, h = 8)

        ont_sum <- go_tree %>%
          group_by(Ontology) %>%
          summarise(n_terms = n_distinct(GOID), .groups = "drop")
        p2 <- ggplot(ont_sum, aes(x = reorder(Ontology, -n_terms),
                                  y = n_terms, fill = Ontology)) +
          geom_col(alpha = 0.85) +
          geom_text(aes(label = comma(n_terms)), vjust = -0.3, size = 3.5) +
          scale_y_continuous(labels = comma) +
          scale_fill_manual(values = c(
            "Biological Process" = "#3574b0",
            "Molecular Function" = "#e8913a",
            "Cellular Component" = "#6a994e")) +
          labs(title = paste0(aid, " — GO ontology composition"),
               x = NULL, y = "Unique GO terms") +
          theme_report + theme(legend.position = "none")
        save_plot(p2, paste0(aid, "_go_ontology_composition"), w = 6, h = 5)
      }
    }
  }

  # ── Top InterPro domains ───────────────────────────────────────────
  ipr_raw <- func$ipr_accessions[func$ipr_accessions != "-" & !is.na(func$ipr_accessions)]
  all_ipr <- trimws(unlist(strsplit(ipr_raw, "[;,]+")))
  all_ipr <- all_ipr[grepl("^IPR", all_ipr)]

  ipr_map <- list()
  for (i in seq_len(nrow(func))) {
    a_str <- func$ipr_accessions[i]; d_str <- func$ipr_descriptions[i]
    if (is.na(a_str) || a_str == "-" || is.na(d_str) || d_str == "-") next
    accs <- trimws(unlist(strsplit(a_str, ";+")))
    descs <- trimws(unlist(strsplit(d_str, ";+")))
    accs <- accs[grepl("^IPR", accs)]
    descs <- descs[descs != "-" & descs != ""]
    for (j in seq_along(accs)) {
      if (!(accs[j] %in% names(ipr_map)) && j <= length(descs)) {
        ipr_map[[accs[j]]] <- descs[j]
      }
    }
  }

  if (length(all_ipr) > 0) {
    ipr_ct <- as.data.frame(table(all_ipr), stringsAsFactors = FALSE)
    names(ipr_ct) <- c("IPR", "Count")
    ipr_ct <- ipr_ct %>% arrange(desc(Count)) %>% slice_head(n = 30)
    ipr_ct$desc <- sapply(ipr_ct$IPR, function(x) ipr_map[[x]] %||% "")
    ipr_ct$label <- ifelse(ipr_ct$desc != "",
      paste0(ipr_ct$IPR, ": ", ifelse(nchar(ipr_ct$desc) > 40,
             paste0(substr(ipr_ct$desc, 1, 37), "..."), ipr_ct$desc)),
      ipr_ct$IPR)

    p <- ggplot(ipr_ct, aes(x = reorder(label, Count), y = Count)) +
      geom_col(fill = "#a855f7", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Top InterPro domains"),
           subtitle = paste0(length(unique(all_ipr)), " unique"),
           x = NULL, y = "Occurrences") + theme_report
    save_plot(p, paste0(aid, "_top_ipr"), w = 11, h = 9)
  }

  # ── Top KEGG terms ─────────────────────────────────────────────────
  kegg_raw <- c(func$eggnog_KEGG_ko, func$kofam_KO)
  kegg_raw <- kegg_raw[kegg_raw != "-" & !is.na(kegg_raw)]
  all_kegg <- trimws(unlist(strsplit(kegg_raw, "[;,]+")))
  all_kegg <- all_kegg[nchar(all_kegg) > 0]

  ko_map <- list()
  for (i in seq_len(nrow(func))) {
    ko <- trimws(func$kofam_KO[i] %||% "-")
    kd <- trimws(func$kofam_def[i] %||% "-")
    if (ko != "-" && kd != "-" && !(ko %in% names(ko_map))) ko_map[[ko]] <- kd
  }

  if (length(all_kegg) > 0) {
    kegg_ct <- as.data.frame(table(all_kegg), stringsAsFactors = FALSE)
    names(kegg_ct) <- c("KO", "Count")
    kegg_ct <- kegg_ct %>% arrange(desc(Count)) %>% slice_head(n = 30)
    kegg_ct$desc <- sapply(kegg_ct$KO, function(x) ko_map[[x]] %||% "")
    kegg_ct$label <- ifelse(kegg_ct$desc != "",
      paste0(kegg_ct$KO, ": ", ifelse(nchar(kegg_ct$desc) > 40,
             paste0(substr(kegg_ct$desc, 1, 37), "..."), kegg_ct$desc)),
      kegg_ct$KO)

    p <- ggplot(kegg_ct, aes(x = reorder(label, Count), y = Count)) +
      geom_col(fill = "#e8913a", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Top KEGG orthologs"),
           subtitle = paste0(length(unique(all_kegg)), " unique"),
           x = NULL, y = "Occurrences") + theme_report
    save_plot(p, paste0(aid, "_top_kegg"), w = 11, h = 9)
  }
}

# ── Per-assembly: top products ───────────────────────────────────────────────
for (gs_path in gene_summary_files) {
  aid <- infer_id(gs_path, c(".gene_product_summary.tsv"))
  gs <- tryCatch(read.delim(gs_path, stringsAsFactors = FALSE),
                 error = function(e) data.frame())
  if (nrow(gs) == 0) next
  top_prod <- gs %>% filter(product != "hypothetical protein") %>%
    count(product, sort = TRUE) %>% slice_head(n = 20)
  if (nrow(top_prod) > 0) {
    p <- ggplot(top_prod, aes(x = reorder(product, n), y = n)) +
      geom_col(fill = "#3574b0", alpha = 0.85) + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Top 20 gene products"),
           x = NULL, y = "Count") + theme_report
    save_plot(p, paste0(aid, "_top_products"), w = 10, h = 7)
  }
}

# ── Per-assembly: telomere density ───────────────────────────────────────────
for (telo_path in telomere_tsv_files) {
  aid <- infer_id(telo_path, c(".telomeres.tsv"))
  telo <- tryCatch(read.delim(telo_path, stringsAsFactors = FALSE),
                   error = function(e) data.frame())
  if (nrow(telo) == 0) next
  names(telo) <- gsub("\\.", "_", tolower(names(telo)))
  if (!"id" %in% names(telo)) names(telo)[1] <- "id"
  fwd <- grep("forward", names(telo), value = TRUE)[1]
  rev_c <- grep("reverse", names(telo), value = TRUE)[1]
  win <- grep("window", names(telo), value = TRUE)[1]
  if (is.na(fwd) || is.na(rev_c) || is.na(win)) next
  telo$total <- telo[[fwd]] + telo[[rev_c]]
  telo$pos_mb <- telo[[win]] / 1e6
  top_s <- telo %>% group_by(id) %>%
    summarise(mx = max(.data[[win]]), .groups = "drop") %>%
    arrange(desc(mx)) %>% slice_head(n = 24) %>% pull(id)
  td <- telo %>% filter(id %in% top_s)
  td$id <- factor(td$id, levels = top_s)
  if (nrow(td) > 0) {
    nr <- ceiling(length(top_s) / 4)
    p <- ggplot(td, aes(x = pos_mb, y = total)) +
      geom_area(fill = "#e05c5c", alpha = 0.5) +
      geom_line(color = "#e05c5c", linewidth = 0.3) +
      facet_wrap(~ id, scales = "free_x", ncol = 4) +
      scale_y_continuous(labels = comma) +
      labs(title = paste0(aid, " — Telomere density"),
           x = "Position (Mb)", y = "Repeats/window") +
      theme_report + theme(strip.text = element_text(size = 7),
                           axis.text = element_text(size = 6))
    save_plot(p, paste0(aid, "_telomere_density"),
              w = 12, h = max(4, nr * 2.2))
  }
}

# ── Per-assembly: telomere positions ─────────────────────────────────────────
for (tg_path in telomere_gff3_files) {
  aid <- infer_id(tg_path, c(".telomeres.gff3"))
  tg <- tryCatch(parse_gff3(tg_path), error = function(e) data.frame())
  if (nrow(tg) == 0) next
  tg$position <- sub(".*position=([^;]+).*", "\\1", tg$attributes)
  tg$position[!grepl("position=", tg$attributes)] <- "unknown"
  tp <- tg %>% count(position)
  p <- ggplot(tp, aes(x = reorder(position, -n), y = n, fill = position)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = comma(n)), vjust = -0.3, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste0(aid, " — Telomere positions"),
         subtitle = paste0(nrow(tg), " regions"),
         x = NULL, y = "Count") +
    theme_report + theme(legend.position = "none")
  p <- smart_y_bar(p, tp$n)
  save_plot(p, paste0(aid, "_telomere_positions"), w = 6, h = 5)
}

# ============================================================================
# ALWAYS-GENERATED PLOTS (BUSCO, OMArk, repeat masking)
# ============================================================================

busco_cols <- c("busco_single","busco_duplicated","busco_fragmented","busco_missing")
if (all(busco_cols %in% names(metrics))) {
  bl <- metrics %>% select(label, all_of(busco_cols)) %>%
    pivot_longer(-label, names_to = "Cat", values_to = "Pct") %>%
    mutate(Cat = factor(recode(Cat,
      busco_single="Single", busco_duplicated="Duplicated",
      busco_fragmented="Fragmented", busco_missing="Missing"),
      levels = rev(c("Single","Duplicated","Fragmented","Missing"))))
  p <- ggplot(bl, aes(x = label, y = Pct, fill = Cat)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = ifelse(Pct >= 3, paste0(round(Pct,1),"%"), "")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(Single="#3574b0", Duplicated="#6ab0e3",
                                 Fragmented="#e8913a", Missing="#d1d5db")) +
    labs(title = "BUSCO completeness", x = NULL, y = "%") + theme_compare
  save_plot(p, "C_busco", w = max(6, 3 + n_distinct(metrics$label) * 1.5))
}

oc <- c("omark_single","omark_dup_unexpected","omark_dup_expected","omark_missing")
if (all(oc %in% names(metrics))) {
  ol <- metrics %>% select(label, all_of(oc)) %>%
    pivot_longer(-label, names_to = "Cat", values_to = "Pct") %>%
    mutate(Cat = factor(recode(Cat,
      omark_single="Single", omark_dup_unexpected="Dup. (unexpected)",
      omark_dup_expected="Dup. (expected)", omark_missing="Missing"),
      levels = rev(c("Single","Dup. (unexpected)","Dup. (expected)","Missing"))))
  p <- ggplot(ol, aes(x = label, y = Pct, fill = Cat)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = ifelse(Pct >= 3, paste0(round(Pct,1),"%"), "")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(Single="#3574b0", `Dup. (unexpected)`="#e8913a",
                                 `Dup. (expected)`="#6ab0e3", Missing="#d1d5db")) +
    labs(title = "OMArk — conserved gene completeness", x = NULL, y = "%") +
    theme_compare
  save_plot(p, "C_omark_completeness",
            w = max(6, 3 + n_distinct(metrics$label) * 1.5))
}

oc2 <- c("omark_consistent","omark_inconsistent","omark_contaminant","omark_unknown")
if (all(oc2 %in% names(metrics))) {
  ol2 <- metrics %>% select(label, all_of(oc2)) %>%
    pivot_longer(-label, names_to = "Cat", values_to = "Pct") %>%
    mutate(Cat = factor(recode(Cat,
      omark_consistent="Consistent", omark_inconsistent="Inconsistent",
      omark_contaminant="Contaminant", omark_unknown="Unknown"),
      levels = rev(c("Consistent","Inconsistent","Contaminant","Unknown"))))
  p <- ggplot(ol2, aes(x = label, y = Pct, fill = Cat)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = ifelse(Pct >= 3, paste0(round(Pct,1),"%"), "")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(Consistent="#6a994e", Inconsistent="#e8913a",
                                 Contaminant="#dc2626", Unknown="#9ca3af")) +
    labs(title = "OMArk — proteome consistency", x = NULL, y = "%") +
    theme_compare
  save_plot(p, "C_omark_consistency",
            w = max(6, 3 + n_distinct(metrics$label) * 1.5))
}

if ("repeat_masked_pct" %in% names(metrics)) {
  p <- ggplot(metrics, aes(x = label, y = repeat_masked_pct)) +
    geom_col(fill = "#8b5e3c", alpha = 0.85) +
    geom_text(aes(label = paste0(round(repeat_masked_pct, 1), "%")),
              vjust = -0.3, size = 3.5) +
    scale_y_continuous(labels = comma) +
    labs(title = "Repeat masking", x = NULL, y = "% genome masked") +
    theme_compare
  save_plot(p, "C_repeat_masking",
            w = max(6, 3 + n_distinct(metrics$label) * 1.5))
}

if ("pct_te_flagged" %in% names(metrics)) {
  p <- ggplot(metrics, aes(x = label, y = pct_te_flagged)) +
    geom_col(fill = "#dc2626", alpha = 0.85) +
    geom_text(aes(label = paste0(round(pct_te_flagged, 1), "%")),
              vjust = -0.3, size = 3.5) +
    scale_y_continuous(labels = comma) +
    labs(title = "TE-flagged genes",
         subtitle = "Putative transposable element-derived gene models",
         x = NULL, y = "% of predicted genes") + theme_compare
  save_plot(p, "C_te_flagged",
            w = max(6, 3 + n_distinct(metrics$label) * 1.5))
}

# ============================================================================
# MULTI-ASSEMBLY ONLY PLOTS
# ============================================================================
if (nrow(metrics) >= 2) {

  feat_cols <- c("n_genes","n_te_flagged","n_trna","n_ncrna","n_telomere","n_centromere")
  feat_cols <- feat_cols[feat_cols %in% names(metrics)]
  if (length(feat_cols) > 0) {
    fl <- metrics %>% select(label, all_of(feat_cols)) %>%
      pivot_longer(-label, names_to = "Feature", values_to = "Count") %>%
      mutate(Feature = recode(Feature, n_genes="Genes", n_te_flagged="TE-flagged",
             n_trna="tRNA", n_ncrna="ncRNA", n_telomere="Telomeres",
             n_centromere="Centromeres"))
    p <- ggplot(fl, aes(x = label, y = Count, fill = Feature)) +
      geom_col(position = "dodge", alpha = 0.85) +
      scale_fill_brewer(palette = "Set2") + scale_y_continuous(labels = comma) +
      labs(title = "Feature counts", x = NULL, y = "Count") + theme_compare
    save_plot(p, "C_feature_counts")
  }

  fc <- c("pct_ipr","pct_eggnog","pct_diamond","pct_kofamscan","pct_deepsig")
  if (all(fc %in% names(metrics))) {
    fl2 <- metrics %>% select(label, all_of(fc)) %>%
      pivot_longer(-label, names_to = "Tool", values_to = "Pct") %>%
      mutate(Tool = recode(Tool, pct_ipr="IPR", pct_eggnog="eggNOG",
             pct_diamond="DIAMOND", pct_kofamscan="KofamScan", pct_deepsig="DeepSig"))
    p <- ggplot(fl2, aes(x = label, y = Pct, fill = Tool)) +
      geom_col(position = "dodge", alpha = 0.85) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Functional tool coverage", x = NULL, y = "% proteins") +
      theme_compare
    save_plot(p, "C_functional_coverage")
  }

  if (all(c("pct_named","pct_hypothetical") %in% names(metrics))) {
    pl <- metrics %>% select(label, pct_named, pct_hypothetical) %>%
      pivot_longer(-label, names_to = "Cat", values_to = "Pct") %>%
      mutate(Cat = recode(Cat, pct_named="Named", pct_hypothetical="Hypothetical"))
    p <- ggplot(pl, aes(x = label, y = Pct, fill = Cat)) +
      geom_col(alpha = 0.85) +
      geom_text(aes(label = paste0(round(Pct,1),"%")),
                position = position_stack(vjust = 0.5), size = 3.5) +
      scale_fill_manual(values = c(Named="#3574b0", Hypothetical="#d1d5db")) +
      labs(title = "Product names", x = NULL, y = "% genes") + theme_compare
    save_plot(p, "C_product_names")
  }

  if (all(c("unique_go_terms","unique_kegg_kos") %in% names(metrics))) {
    gl <- metrics %>% select(label, unique_go_terms, unique_kegg_kos) %>%
      pivot_longer(-label, names_to = "Metric", values_to = "Count") %>%
      mutate(Metric = recode(Metric,
             unique_go_terms="GO terms", unique_kegg_kos="KEGG KOs"))
    p <- ggplot(gl, aes(x = label, y = Count, fill = Metric)) +
      geom_col(position = "dodge", alpha = 0.85) +
      scale_fill_manual(values = c("GO terms"="#6a994e","KEGG KOs"="#e8913a")) +
      scale_y_continuous(labels = comma) +
      labs(title = "Functional vocabulary", x = NULL, y = "Count") + theme_compare
    save_plot(p, "C_go_kegg_diversity")
  }

  rc <- c("busco_complete","omark_consistent","pct_named",
          "repeat_masked_pct","pct_diamond","pct_te_flagged")
  rc <- rc[rc %in% names(metrics)]
  if (length(rc) >= 3) {
    # Invert TE-flagged so lower = better displays correctly
    rl <- metrics %>% select(label, all_of(rc)) %>%
      mutate(across(any_of("pct_te_flagged"), ~ 100 - .x)) %>%
      pivot_longer(-label, names_to = "Metric", values_to = "Score") %>%
      mutate(Metric = recode(Metric, busco_complete="BUSCO",
             omark_consistent="OMArk", pct_named="Named",
             repeat_masked_pct="Masked", pct_diamond="SwissProt",
             pct_te_flagged="Non-TE (100-TE%)"))
    p <- ggplot(rl, aes(x = Metric, y = Score, color = label, group = label)) +
      geom_point(size = 3) + geom_line(linewidth = 0.8, alpha = 0.7) +
      scale_color_brewer(palette = "Set1") + coord_flip() +
      labs(title = "Quality summary", x = NULL, y = "Score (%)",
           color = "Assembly") + theme_compare
    save_plot(p, "C_quality_summary")
  }
}

# ── Done ─────────────────────────────────────────────────────────────────────
n_plots <- length(list.files(opt$outdir, pattern = "\\.png$"))
message(sprintf("Generated %d plots in %s/", n_plots, opt$outdir))