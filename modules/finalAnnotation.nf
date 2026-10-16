process FINALIZE_ANNOTATION {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/final", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(mode), path(original_files, stageAs: 'original/*'), path(validation), path(gff, stageAs: 'refined.gff3'), path(utr_status)
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: bundle
    script:
    """
    python final_annotation_inputs.py --assembly-id '${aid}' --mode '${mode}' --prediction original \
      --validation '${validation}' --gff '${gff}' --utr-status '${utr_status}' --outdir '${aid}'
    """
}

process REPORT_FINAL_AED {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/aed/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(gff), path(identity), path(indexed), path(annotated)
    path helpers
    output:
    tuple val(aid), path('aed_scores.tsv'), emit: scores
    tuple val(aid), path('aed_status.json'), emit: status
    tuple val(aid), path('README.md'), emit: report
    script:
    """
    python final_annotation.py aed --gff '${gff}' --identity '${identity}' \
      --inputs '${indexed}/inputs.json' --annotated '${annotated}' \
      --cds-only '${params.ingenannot_rna_cds_only}' --penalty ${params.ingenannot_rna_penalty} --outdir .
    """
}
