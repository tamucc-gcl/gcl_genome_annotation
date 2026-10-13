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
    publishDir "${params.outdir}/qc/baseline/${assembly_id}/omark", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), val(taxid), path(protein_fasta), path(omamer_db), path(taxa_db)
    
    output:
    tuple val(assembly_id), path("${assembly_id}.png"), emit: omark_plot
    tuple val(assembly_id), path("${assembly_id}_detailed.txt"), emit: omark_detailed
    tuple val(assembly_id), path("${assembly_id}.summary.txt"), emit: omark_summary
    tuple val(assembly_id), path("${assembly_id}.ump"),          emit: omark_ump
    tuple val(assembly_id), path('omark_output'), emit: results
    
    script:
    """
    set -euo pipefail
    omamer search \
        --nthreads ${task.cpus} \\
        --db ${omamer_db} \\
        --query ${protein_fasta} \\
        --out ${assembly_id}.omamer

    omark \
        --file ${assembly_id}.omamer \\
        --database ${omamer_db} \\
        --ete_ncbi_db ${taxa_db} \\
        --taxid '${taxid}' \\
        --og_fasta ${protein_fasta}
    
    cp omark_output/${assembly_id}.png ${assembly_id}.png
    cp omark_output/${assembly_id}_detailed_summary.txt ${assembly_id}_detailed.txt
    cp omark_output/${assembly_id}.sum ${assembly_id}.summary.txt
    cp omark_output/${assembly_id}.ump ${assembly_id}.ump
    # OMArk's version flag still requires analysis inputs in some releases.
    # Read package metadata without invoking its argument parser.
    python -c "from importlib.metadata import version; print(version('omark'))" > omark_output/version.txt
    python -c "from importlib.metadata import version; print(version('omamer'))" > omark_output/omamer_version.txt
    """
    
    stub:
    """
    touch ${assembly_id}.png
    touch ${assembly_id}_detailed.txt
    touch ${assembly_id}.summary.txt
    touch ${assembly_id}.ump
    mkdir -p omark_output
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
    tuple val(oma_db_name), path("${oma_db_name}.h5"), path('taxa.sqlite'), emit: databases
    
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
