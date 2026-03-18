/*
========================================================================================
    ANNOTATION REPORT MODULE
========================================================================================
    Generate a comprehensive Markdown annotation report with embedded plots.

    Collects outputs from all upstream modules and produces:
      - PNG plots via R script (annotation_report_plots.R)
      - Markdown report via Python script (generate_annotation_report.py)
      - A plots/ directory with all generated figures

    The report is a standalone Markdown file viewable on GitHub, with
    relative image references to the plots/ directory.

    Inputs:
        - assembly_id:          Sample identifier
        - assembly_stats:       Assembly stats TSV (from FILTER_ASSEMBLY)
        - annotated_gff3:       Final merged annotated GFF3
        - gene_summary:         Gene product summary TSV
        - functional_tsv:       Merged functional annotations TSV
        - repeat_tbl:           RepeatMasker summary table
        - repeat_tbl_full:      buildSummary.pl detailed table
        - repeat_landscape:     Repeat landscape HTML
        - repeat_divsum:        Divergence summary file
        - busco_txt:            BUSCO short summary
        - omark_summary:        OMArk summary
        - omark_detailed:       OMArk detailed summary
        - agat_stats:           AGAT annotation statistics
        - trnascan_tsv:         tRNAscan-SE tabular output
        - infernal_tblout:      Infernal tblout (merged)
        - telomere_gff3:        Telomere annotations GFF3
        - telomere_svg:         Telomere density plot (tidk)
        - centromere_gff3:      Centromere annotations GFF3

    Outputs:
        - report_md:    ${assembly_id}_annotation_report.md
        - plots_dir:    ${assembly_id}_annotation_plots/

    Tools: R (ggplot2, dplyr, tidyr, scales, patchwork, treemapify),
           Python 3 (stdlib only)
========================================================================================
*/

process ANNOTATION_REPORT {
    tag "${assembly_id}"
    label 'annotation_report'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id),
          path(assembly_stats),
          path(annotated_gff3),
          path(gene_summary),
          path(functional_tsv),
          path(repeat_tbl),
          path(repeat_tbl_full),
          path(repeat_landscape),
          path(repeat_divsum),
          path(busco_txt),
          path(omark_summary),
          path(omark_detailed),
          path(agat_stats),
          path(trnascan_tsv),
          path(infernal_tblout),
          path(telomere_gff3),
          path(telomere_svg),
          path(centromere_gff3)

    output:
    tuple val(assembly_id),
          path("${assembly_id}_annotation_report.md"),
          emit: report_md
    tuple val(assembly_id),
          path("${assembly_id}_annotation_plots"),
          emit: plots_dir

    script:
    """
    set -euo pipefail

    mkdir -p ${assembly_id}_annotation_plots

    # ── 1. Generate plots with R ──────────────────────────────────────────────
    Rscript ${projectDir}/r_scripts/annotation_report_plots.R \
        --assembly_id ${assembly_id} \
        --gff3 ${annotated_gff3} \
        --gene_summary ${gene_summary} \
        --functional_tsv ${functional_tsv} \
        --repeat_tbl ${repeat_tbl} \
        --repeat_tbl_full ${repeat_tbl_full} \
        --repeat_divsum ${repeat_divsum} \
        --busco_txt ${busco_txt} \
        --trnascan_tsv ${trnascan_tsv} \
        --infernal_tblout ${infernal_tblout} \
        --telomere_gff3 ${telomere_gff3} \
        --centromere_gff3 ${centromere_gff3} \
        --outdir ${assembly_id}_annotation_plots

    # ── 2. Assemble Markdown report ───────────────────────────────────────────
    python3 ${projectDir}/py_scripts/generate_annotation_report.py \
        --assembly_id ${assembly_id} \
        --assembly_stats ${assembly_stats} \
        --gff3 ${annotated_gff3} \
        --gene_summary ${gene_summary} \
        --functional_tsv ${functional_tsv} \
        --repeat_tbl ${repeat_tbl} \
        --repeat_tbl_full ${repeat_tbl_full} \
        --busco_txt ${busco_txt} \
        --omark_summary ${omark_summary} \
        --omark_detailed ${omark_detailed} \
        --agat_stats ${agat_stats} \
        --trnascan_tsv ${trnascan_tsv} \
        --infernal_tblout ${infernal_tblout} \
        --telomere_gff3 ${telomere_gff3} \
        --centromere_gff3 ${centromere_gff3} \
        --plots_dir ${assembly_id}_annotation_plots \
        --outfile ${assembly_id}_annotation_report.md
    """

    stub:
    """
    mkdir -p ${assembly_id}_annotation_plots
    touch ${assembly_id}_annotation_report.md
    """
}