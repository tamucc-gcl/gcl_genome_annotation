/*
========================================================================================
    OMArk MODULE
========================================================================================
    Assess genome completeness using OMArk
========================================================================================
*/

process OMARK {
    tag "${assembly_id}"
    label 'omark'
    publishDir "${params.outdir}/qc/omark", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), path(protein_fasta)
    path(omamer_db)
    path(taxa_db)
    
    output:
    tuple val(assembly_id), path("${assembly_id}.png"), emit: omark_plot
    tuple val(assembly_id), path("${assembly_id}_detailed.txt"), emit: omark_detailed
    tuple val(assembly_id), path("${assembly_id}.summary.txt"), emit: omark_summary
    
    script:
    """
    omamer search \
        --nthreads ${task.cpus} \\
        --db ${omamer_db} \\
        --query ${protein_fasta} \\
        --out ${assembly_id}.omamer

    omark \
        --file ${assembly_id}.omamer \\
        --database ${omamer_db} \\
        --ete_ncbi_db ${taxa_db} \\
        --taxid ${params.ncbi_speciesID} \\
        --og_fasta ${protein_fasta}
    
    mv omark_output/${assembly_id}.png ${assembly_id}.png
    mv omark_output/${assembly_id}_detailed_summary.txt ${assembly_id}_detailed.txt
    mv omark_output/${assembly_id}.sum ${assembly_id}.summary.txt
    """
    
    stub:
    """
    touch ${assembly_id}.png
    touch ${assembly_id}_detailed_summary.txt
    touch ${assembly_id}.sum
    """
}

/*
========================================================================================
    DOWNLOAD OMArk DATABASES MODULE
========================================================================================
    Downloads the specified databases needed for OMArk once, before parallel OMArk runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_OMARK_DB {
    tag "${oma_db_name}"
    label 'omark_download'
    
    storeDir "${params.omark_database}"
    
    input:
    val oma_db_name

    output:
    path("${oma_db_name}.h5"), emit: omamer_db
    path("taxa.sqlite"), emit: taxa_db
    
    script:
    """
    #Download oma LUCA database
    wget https://omabrowser.org/All/${oma_db_name}.h5

    #Download proper format NCBI taxonomy
    python - <<PY
    from ete3 import NCBITaxa
    dbfile = "./taxa.sqlite"
    ncbi = NCBITaxa(dbfile=dbfile)
    print("dbfile BEFORE:", ncbi.dbfile)
    ncbi.update_taxonomy_database()
    print("dbfile AFTER:", ncbi.dbfile)
    PY

    """
    
    stub:
    """
    touch ${params.oma_db_name}.h5
    touch taxa.sqlite
    """
}