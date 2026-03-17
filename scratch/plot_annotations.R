library(Biostrings)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)

#---------------------------
# input files
#---------------------------
fasta_file <- "~/../Downloads/Sde-CMat_203_hap2.fasta"
telomere_gff_file   <- "~/../Downloads/Sde-CMat_203_hap2.telomeres.gff3"
centromere_gff_file <- '~/../Downloads/Sde-CMat_203_hap2.centromere.gff3'
repeat_gff_file <- '~/../Downloads/Sde-CMat_203_hap2.gff'

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
                type) %>%
  mutate(width = end - start + 1)

centromere_gff <- import(centromere_gff_file) %>%
  as_tibble() %>%
  dplyr::select(seqid = seqnames,
                start,
                end,
                strand,
                type) %>%
  mutate(width = end - start + 1)

repeat_gff <- import(repeat_gff_file) %>%
  as_tibble() %>%
  mutate(repeat_name = str_extract(Target, "^[^ ]+"))
  
repeat_df <- repeat_gff %>%
  mutate(repeat_class = case_when(
    str_detect(repeat_name, regex("^ltr|_ltr|gypsy|copia|erv", ignore_case = TRUE)) ~ "LTR",
    str_detect(repeat_name, regex("line|l1|l2|cr1|rex", ignore_case = TRUE)) ~ "LINE",
    str_detect(repeat_name, regex("sine|alu|mir", ignore_case = TRUE)) ~ "SINE",
    str_detect(repeat_name, regex("dna|hat|helitron|mariner|tc1|piggybac|mule", ignore_case = TRUE)) ~ "DNA_transposon",
    str_detect(repeat_name, regex("sat|satellite|cent", ignore_case = TRUE)) ~ "Satellite",
    str_detect(repeat_name, regex("simple|microsatellite|tandem", ignore_case = TRUE)) ~ "Simple_repeat",
    str_detect(repeat_name, regex("low_complexity", ignore_case = TRUE)) ~ "Low_complexity",
    str_detect(repeat_name, regex("rrna|trna|snrna|srp", ignore_case = TRUE)) ~ "RNA_repeat",
    str_detect(repeat_name, regex("^rnd", ignore_case = TRUE)) ~ "Unknown_de_novo",
    str_detect(repeat_name, regex("^DR[0-9]+$", ignore_case = TRUE)) ~ "Repbase_or_known",
    TRUE ~ "Other"))

  dplyr::select(seqid = seqnames,
                start,
                end,
                strand,
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
  geom_rect(data = filter(plot_gff, type != 'dispersed_repeat'),
            aes(xmin = start,
                xmax = end,
                ymin = as.numeric(seqid) - 0.35,
                ymax = as.numeric(seqid) + 0.35,
                fill = type),
            color = NA) +
  # scale_y_continuous(breaks = levels(plot_scaffolds$seqid),
  #                    labels = levels(plot_scaffolds$seqid)) +
  scale_x_continuous(labels = scales::comma_format(scale = 1/1e6, suffix = 'M')) +
  labs(x = "Position on scaffold",
       y = "Scaffold",
       fill = "Annotation type") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

tel <- gff_df %>%
  filter(type == "telomeric_repeat") %>%
  left_join(scaffold_df, by = "seqid") %>%
  mutate(
    dist_to_start = start,
    dist_to_end = length - end,
    near_end = dist_to_start < 50000 | dist_to_end < 50000
  )

tel %>%
  group_by(near_end) %>%
  summarise(
    n = n(),
    mean_width = mean(end - start)
  )

tel %>%
  mutate(width = end - start) %>%
  ggplot(aes(x = pmin(dist_to_start, dist_to_end), y = width)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Distance to nearest scaffold end", y = "Telomere length")
