// Existing reference entry point, extended with per-assembly settings and integrity artifacts.
process FILTER_ASSEMBLY {
    tag "${assembly_id}"
    label 'filter_assembly'
    publishDir "${params.outdir}/filtered_assemblies", mode: params.publish_dir_mode
    input:
    tuple val(assembly_id), path(assembly_path), val(min_contig_length)
    path helper
    output:
    tuple val(assembly_id), path("${assembly_id}"), emit: bundle
    script:
    """
    python ${helper} reference --source '${assembly_path}' --assembly ${assembly_id} --minimum ${min_contig_length} --outdir ${assembly_id}
    """
}
