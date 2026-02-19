/*
========================================================================================
    BUSCO MODULE
========================================================================================
    Assess genome completeness using BUSCO
    Note: Set params.busco_lineage in config (e.g., "actinopterygii_odb10")
========================================================================================
*/

process BUSCO {
    tag "${assembly_id}"
    label 'busco'
    publishDir "${params.outdir}/qc/busco", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), path(protein_fasta)
    path(busco_db)
    
    output:
    tuple val(assembly_id), path("${assembly_id}_busco.txt"), emit: busco
    
    script:
    def lineage = params.busco_lineage ?: "auto"
    """
    busco \\
        --in ${protein_fasta} \\
        --out ${assembly_id}_busco \\
        --mode proteins \\
        --lineage_dataset ${lineage} \\
        --cpu ${task.cpus} \\
        --offline \\
        --download_path ${params.busco_downloads ?: '/path/to/busco_downloads'}
    mv ${assembly_id}_busco/short_summary.specific.${lineage}.${assembly_id}_busco.txt ${assembly_id}_busco.txt 
    """
    
    stub:
    """
    touch ${assembly_id}_busco.txt
    """
}