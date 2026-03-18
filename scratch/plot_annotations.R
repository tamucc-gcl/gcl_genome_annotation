library(Biostrings)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(legendry)
library(purrr)

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
gff_df <- map_dfr(c(telomere_gff_file,
                    centromere_gff_file,
                    repeat_gff_file,
                    ncrna_gff_file),
                  
                  ~import(.x) %>%
                    as_tibble() %>%
                    dplyr::select(seqid = seqnames,
                                  start,
                                  end,
                                  strand,
                                  source,
                                  type) %>%
                    mutate(width = end - start + 1))

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
  geom_rect(data = filter(plot_gff, width > 1000 | 
                            source %in% c('tidk', 'centroAnno', 'tRNAscan-SE')),
            aes(xmin = start,
                xmax = end,
                ymin = as.numeric(seqid) - 0.35,
                ymax = as.numeric(seqid) + 0.35,
                fill = interaction(source, type)),
            color = NA, alpha = 1) +
  # scale_y_continuous(breaks = levels(plot_scaffolds$seqid),
  #                    labels = levels(plot_scaffolds$seqid)) +
  scale_x_continuous(labels = scales::comma_format(scale = 1/1e6, suffix = 'M')) +
  facet_wrap( ~ source) +
  guides(fill = guide_legend_group()) +
  labs(x = "Position on scaffold",
       y = "Scaffold",
       fill = "Annotation type") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

