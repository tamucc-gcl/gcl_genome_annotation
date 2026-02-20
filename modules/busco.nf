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

/*
========================================================================================
    DOWNLOAD BUSCO DATABASE MODULE
========================================================================================
    Downloads the specified BUSCO lineage database once, before parallel BUSCO runs.
    This prevents multiple simultaneous download attempts when running BUSCO in parallel.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_BUSCO_DB {
    tag "${lineage}"
    label 'busco_download'
    
    storeDir "${params.busco_downloads}"
    
    input:
    val(lineage)
    
    output:
    path("lineages/${lineage}"), emit: db
    
    script:
    """
    # Create lineages subdirectory structure expected by BUSCO
    mkdir -p lineages
    
    # Download the specified lineage dataset
    busco --download ${lineage} --download_path .
    
    # Verify download succeeded
    if [ ! -d "lineages/${lineage}" ]; then
        echo "ERROR: Failed to download BUSCO lineage: ${lineage}"
        exit 1
    fi
    
    echo "Successfully downloaded BUSCO lineage: ${lineage}"
    """
    
    stub:
    """
    mkdir -p lineages/${lineage}
    touch lineages/${lineage}/dataset.cfg
    """
}