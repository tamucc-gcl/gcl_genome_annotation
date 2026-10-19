// Existing FILTER_GENE_MODELS, adapted to the accepted assessment handoff.
// Native Nextflow staging/cache; selection writes distinct retained/excluded sets.
process FILTER_GENE_MODELS {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/selection", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(gff), path(proteins), path(cds), path(validation), path(identity), path(assessment), path(assessment_status)
    val selection_settings
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: bundle
    script:
    """
    python model_selection.py --assembly-id '${aid}' --gff '${gff}' --proteins '${proteins}' --cds '${cds}' \
      --validation '${validation}' --identity '${identity}' --assessment '${assessment}' \
      --assessment-status '${assessment_status}' --settings '${selection_settings}' --outdir '${aid}'
    """
}
