/*
========================================================================================
    ANNOTATION REPORT MODULE
========================================================================================
    Generate a unified Markdown annotation report covering ALL assemblies,
    with per-assembly sections and cross-assembly comparisons.

    This process runs ONCE after all per-assembly annotations complete.
    All inputs are collected (list) channels — one file per assembly.
    Assembly IDs are inferred from filenames.

    Inputs (all collected paths — one file per assembly):
        - assembly_stats:    Assembly stats TSVs
        - annotated_gff3s:   Annotated GFF3 files
        - gene_summaries:    Gene product summary TSVs
        - functional_tsvs:   Merged functional annotation TSVs
        - repeat_tbls:       RepeatMasker summary tables
        - repeat_tbls_full:  buildSummary.pl detailed tables
        - busco_txts:        BUSCO short summaries
        - omark_summaries:   OMArk summary files
        - omark_detaileds:   OMArk detailed files
        - agat_stats:        AGAT annotation statistics
        - trnascan_tsvs:     tRNAscan-SE tabular outputs
        - infernal_tblouts:  Infernal tblouts (merged)
        - telomere_gff3s:    Telomere annotation GFF3s
        - telomere_tsvs:     tidk windowed telomere count TSVs
        - centromere_gff3s:   Centromere annotation GFF3s

    Outputs:
        - report_md:    annotation_report.md
        - plots_dir:    annotation_report_plots/
        - metrics_tsv:  annotation_metrics.tsv

    Tools: R (ggplot2, dplyr, tidyr, scales, patchwork, optparse),
           Python 3 (stdlib only)
========================================================================================
*/

process ANNOTATION_REPORT {
    tag "report"
    label 'annotation_report'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    input:
    path(assembly_stats)
    path(annotated_gff3s)
    path(gene_summaries)
    path(functional_tsvs)
    path(repeat_tbls)
    path(busco_txts)
    path(omark_summaries)
    path(omark_detaileds)
    path(agat_stats)
    path(trnascan_tsvs)
    path(infernal_tblouts)
    path(telomere_gff3s)
    path(telomere_tsvs)
    path(centromere_gff3s)

    output:
    path("annotation_report.md"),       emit: report_md
    path("annotation_report_plots"),    emit: plots_dir
    path("annotation_metrics.tsv"),     emit: metrics_tsv

    script:
    """
    set -euo pipefail

    mkdir -p annotation_report_plots

    # ── 1. Collect per-assembly metrics into a single TSV ─────────────────────
    python3 ${projectDir}/scripts/collect_annotation_metrics.py \
        --gff3s ${annotated_gff3s} \
        --gene_summaries ${gene_summaries} \
        --functional_tsvs ${functional_tsvs} \
        --repeat_tbls ${repeat_tbls} \
        --busco_txts ${busco_txts} \
        --omark_summaries ${omark_summaries} \
        --outfile annotation_metrics.tsv

    # ── 2. Generate plots with R ──────────────────────────────────────────────
    Rscript ${projectDir}/scripts/annotation_report_plots.R \
        --metrics annotation_metrics.tsv \
        --gff3s ${annotated_gff3s} \
        --gene_summaries ${gene_summaries} \
        --functional_tsvs ${functional_tsvs} \
        --busco_txts ${busco_txts} \
        --omark_summaries ${omark_summaries} \
        --trnascan_tsvs ${trnascan_tsvs} \
        --infernal_tblouts ${infernal_tblouts} \
        --telomere_gff3s ${telomere_gff3s} \
        --telomere_tsvs ${telomere_tsvs} \
        --centromere_gff3s ${centromere_gff3s} \
        --outdir annotation_report_plots

    # ── 3. Assemble Markdown report ───────────────────────────────────────────
    python3 ${projectDir}/scripts/generate_annotation_report.py \
        --metrics annotation_metrics.tsv \
        --assembly_stats ${assembly_stats} \
        --repeat_tbls ${repeat_tbls} \
        --busco_txts ${busco_txts} \
        --omark_summaries ${omark_summaries} \
        --omark_detaileds ${omark_detaileds} \
        --agat_stats ${agat_stats} \
        --plots_dir annotation_report_plots \
        --outfile annotation_report.md
    """

    stub:
    """
    mkdir -p annotation_report_plots
    touch annotation_report.md
    touch annotation_metrics.tsv
    """
}