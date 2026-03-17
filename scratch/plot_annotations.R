library(Biostrings)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(legendry)

#---------------------------
# input files
#---------------------------
fasta_file <- "~/../Downloads/Sde-CMat_203_hap2.fasta"
telomere_gff_file   <- "~/../Downloads/Sde-CMat_203_hap2.telomeres.gff3"
centromere_gff_file <- '~/../Downloads/Sde-CMat_203_hap2.centromere.gff3'
repeat_gff_file <- '~/../Downloads/Sde-CMat_203_hap2.gff'
ncrna_gff_file <- '~/../Downloads/Sde-CMat_203_hap2.trnascan.gff3'

#---------------------------
# read fasta and get scaffold lengths
#---------------------------
fa <- readDNAStringSet(fasta_file)

scaffold_df <- tibble(seqid = names(fa),
                      length = width(fa))


#---------------------------
# read gff3(s)
#---------------------------
telomere_gff <- import(telomere_gff_file) %>%
  as_tibble() %>%
  dplyr::select(seqid = seqnames,
                start,
                end,
                strand,
                source,
                type) %>%
  mutate(width = end - start + 1)

centromere_gff <- import(centromere_gff_file) %>%
  as_tibble() %>%
  dplyr::select(seqid = seqnames,
                start,
                end,
                strand,
                source,
                type) %>%
  mutate(width = end - start + 1)

repeat_gff <- import(repeat_gff_file) %>%
  as_tibble() %>%
  dplyr::select(seqid = seqnames,
                start,
                end,
                strand,
                source,
                type) %>%
  mutate(width = end - start + 1)

gff_df <- bind_rows(telomere_gff,
                    centromere_gff,
                    repeat_gff)

#---------------------------
# keep only top N scaffolds that have annotations
#---------------------------
plot_scaffolds <- scaffold_df %>%
  semi_join(gff_df, by = "seqid") %>%
  slice_max(length, n = 16) %>%
  mutate(seqid = fct_reorder(seqid, length))

plot_gff <- gff_df %>%
  semi_join(plot_scaffolds, by = "seqid") %>%
  mutate(seqid = factor(seqid, levels = levels(plot_scaffolds$seqid)))

#---------------------------
# plot
#---------------------------
ggplot() +
  # scaffold background
  geom_segment(data = plot_scaffolds,
               aes(x = 1, 
                   xend = length, 
                   y = seqid, 
                   yend = seqid),
               linewidth = 2,
               color = "gray75") +
  # annotation blocks
  geom_rect(data = filter(plot_gff, width > 1000 | source %in% c('tidk', 'centroAnno')),
            aes(xmin = start,
                xmax = end,
                ymin = as.numeric(seqid) - 0.35,
                ymax = as.numeric(seqid) + 0.35,
                fill = interaction(source, type)),
            color = NA, alpha = 0.5) +
  # scale_y_continuous(breaks = levels(plot_scaffolds$seqid),
  #                    labels = levels(plot_scaffolds$seqid)) +
  scale_x_continuous(labels = scales::comma_format(scale = 1/1e6, suffix = 'M')) +
  guides(fill = guide_legend_group()) +
  labs(x = "Position on scaffold",
       y = "Scaffold",
       fill = "Annotation type") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

