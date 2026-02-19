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
    tag "${lineage}"
    label 'omark_download'
    
    storeDir "${params.omark_downloads}"
    
    input:
    
    output:
    path("${params.oma_db_name}.h5"), emit: omamer_db
    path("taxa.sqlite"), emit: taxa_db
    
    script:
    """
    #Download oma LUCA database
    wget https://omabrowser.org/All/${params.oma_db_name}.h5

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