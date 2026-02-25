/*
========================================================================================
    Infernal MODULE
========================================================================================
    Run Infernal on masked genome
========================================================================================
*/

process INFERNAL {
    tag "${assembly_id}"
    label 'infernal'
    publishDir "${params.outdir}/ncRNA", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(masked_assembly)
    path rfam_database

    output:
    tuple val(assembly_id), path("${assembly_id}_rfam.tblout"),            emit: infernal_tblout
    tuple val(assembly_id), path("${assembly_id}_rfam.cmscan"),            emit: infernal_cmscan

    script:
    """
    cmscan \
        --cpu ${task.cpus} \
        --cut_ga \
        --rfam \
        --nohmmonly \
        --tblout ${assembly_id}_rfam.tblout \
        --fmt 2 \
        --clanin ${rfam_database}/Rfam.clanin \
        ${rfam_database}/Rfam.cm \
        ${masked_assembly} > ${assembly_id}_rfam.cmscan
    """

    stub:
    """
    touch ${assembly_id}_rfam.tblout
    touch ${assembly_id}_rfam.cmscan
    """
}
/*
========================================================================================
    DOWNLOAD RFam DATABASE MODULE
========================================================================================
    Downloads the specified Rfam database once, before parallel Rfam runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_RFAM_DB {
    tag "${rfam_version}"
    label 'download_rfam_db'
    storeDir "${params.rfam_database}"

    input:
    val rfam_version

    output:
    path("rfam_db"), emit: rfam_db_dir

    script:
    """
    set -euo pipefail

    mkdir rfam_db
    cd rfam_db
    wget --tries=10 --timeout=30 https://ftp.ebi.ac.uk/pub/databases/Rfam/${rfam_version}/Rfam.cm.gz
    wget --tries=10 --timeout=30 https://ftp.ebi.ac.uk/pub/databases/Rfam/${rfam_version}/Rfam.clanin
    gunzip Rfam.cm.gz
    cmpress Rfam.cm
    """
}