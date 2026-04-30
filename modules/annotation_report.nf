/*
========================================================================================
    ANNOTATION REPORT MODULE
========================================================================================
    Generate a unified Markdown annotation report covering ALL assemblies.
    Runs ONCE after all per-assembly annotations complete.
    All inputs are collected (list) channels — one file per assembly.
    Assembly IDs are inferred from filenames.

    Inputs (all collected paths — one file per assembly):
        - assembly_stats:       Assembly stats TSVs
        - annotated_gff3s:      Annotated GFF3 files
        - gene_summaries:       Gene product summary TSVs
        - functional_tsvs:      Merged functional annotation TSVs
        - repeat_tbls:          RepeatMasker summary tables
        - repeat_tbls_full:     buildSummary.pl detailed tables
        - busco_txts:           BUSCO short summaries
        - omark_summaries:      OMArk summary files
        - omark_detaileds:      OMArk detailed files
        - agat_stats:           AGAT annotation statistics
        - trnascan_tsvs:        tRNAscan-SE tabular outputs
        - infernal_tblouts:     Infernal tblouts (merged)
        - telomere_gff3s:       Telomere annotation GFF3s
        - telomere_tsvs:        tidk windowed telomere count TSVs
        - centromere_gff3s:     Centromere annotation GFF3s
        - structural_qc_tsvs:   Per-assembly structural QC summary TSVs

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
    publishDir "${params.outdir}",         mode: params.publish_dir_mode, pattern: "annotation_report.md"
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode, pattern: "annotation_report_plots"
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode, pattern: "annotation_metrics.tsv"

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
    path(structural_qc_tsvs)

    output:
    path("annotation_report.md"),       emit: report_md
    path("annotation_report_plots"),    emit: plots_dir
    path("annotation_metrics.tsv"),     emit: metrics_tsv

    script:
    // Serialise all relevant params into a JSON-like string for the report
    def db_info = [
        dfam_version:             params.dfam_version             ?: "unset",
        taxa_name:                params.taxa_name                ?: "unset",
        orthodb_clade:            params.orthodb_clade            ?: "none",
        ref_proteome_accession:   params.ref_proteome_accession   ?: "none",
        busco_lineage:            params.busco_lineage            ?: "unset",
        oma_db_name:              params.oma_db_name              ?: "unset",
        rfam_version:             params.rfam_version             ?: "unset",
        interproscan_version:     params.interproscan_version     ?: "unset",
        eggnog_version:           params.eggnog_version           ?: "unset",
        diamond_version:          params.diamond_version          ?: "unset",
        kofamscan_version:        params.kofamscan_version        ?: "unset",
        telomere_motif:           params.telomere_motif           ?: "TTAGGG",
        ncbi_speciesID:           params.ncbi_speciesID           ?: "none",
        braker_rounds:            params.braker_rounds            ?: 8,
        n_rm_chunks:              params.n_rm_chunks              ?: 8,
        filter_min_cds_length:    params.filter_min_cds_length    ?: 300,
        filter_min_gene_length:   params.filter_min_gene_length   ?: 200,
        filter_min_intron_length: params.filter_min_intron_length ?: 10,
        te_flag_overlap_frac:     params.te_flag_overlap_frac     ?: 0.5,
        min_contig_length:        params.min_contig_length        ?: 50000,
        terrier_threshold:        params.terrier_threshold        ?: 0.7,
        repeat_filter_pident:     params.repeat_filter_pident     ?: 25,
    ]
    // Build shell variable as TSV pairs consumed by generate_annotation_report.py
    def db_info_str = db_info.collect { k, v -> "${k}\\t${v}" }.join("\\n")
    """
    set -euo pipefail

    mkdir -p annotation_report_plots

    # Write database/parameter info file for report script
    printf '${db_info_str}\\n' > pipeline_params.tsv

    # ── 1. Collect per-assembly metrics into a single TSV ─────────────────────
    python3 ${projectDir}/py_scripts/collect_annotation_metrics.py \
        --gff3s ${annotated_gff3s} \
        --gene_summaries ${gene_summaries} \
        --functional_tsvs ${functional_tsvs} \
        --repeat_tbls ${repeat_tbls} \
        --busco_txts ${busco_txts} \
        --omark_summaries ${omark_summaries} \
        --structural_qc_tsvs ${structural_qc_tsvs} \
        --outfile annotation_metrics.tsv

    # ── 2. Generate plots with R ──────────────────────────────────────────────
    Rscript ${projectDir}/r_scripts/annotation_report_plots.R \
        --metrics annotation_metrics.tsv \
        --gff3s "${annotated_gff3s}" \
        --gene_summaries "${gene_summaries}" \
        --functional_tsvs "${functional_tsvs}" \
        --busco_txts "${busco_txts}" \
        --omark_summaries "${omark_summaries}" \
        --trnascan_tsvs "${trnascan_tsvs}" \
        --infernal_tblouts "${infernal_tblouts}" \
        --telomere_gff3s "${telomere_gff3s}" \
        --telomere_tsvs "${telomere_tsvs}" \
        --centromere_gff3s "${centromere_gff3s}" \
        --outdir annotation_report_plots

    # ── 3. Assemble Markdown report ───────────────────────────────────────────
    python3 ${projectDir}/py_scripts/generate_annotation_report.py \
        --metrics annotation_metrics.tsv \
        --assembly_stats ${assembly_stats} \
        --repeat_tbls ${repeat_tbls} \
        --busco_txts ${busco_txts} \
        --omark_summaries ${omark_summaries} \
        --omark_detaileds ${omark_detaileds} \
        --agat_stats ${agat_stats} \
        --structural_qc_tsvs ${structural_qc_tsvs} \
        --pipeline_params pipeline_params.tsv \
        --plots_dir annotation_report_plots \
        --plots_md_path reports/annotation_report_plots \
        --outfile annotation_report.md
    """

    stub:
    """
    mkdir -p annotation_report_plots
    touch annotation_report.md
    touch annotation_metrics.tsv
    """
}