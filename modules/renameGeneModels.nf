// Adapt the existing assembly_gN/.tN scheme to validated selected models.
process RENAME_GENE_MODELS {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/identified", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(selection,stageAs:'selected_source'), path(identity)
    val id_prefix
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: bundle
    script:
    def prefixOption = id_prefix ? "--prefix '${id_prefix}_${aid}'" : ''
    """
    python final_identifiers.py --assembly-id '${aid}' --selection '${selection}' --identity '${identity}' ${prefixOption} --outdir '${aid}'
    """
}
