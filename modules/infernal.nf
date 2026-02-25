/*
========================================================================================
    Infernal MODULE
========================================================================================
    Run Infernal cmscan on a genome chunk (or full genome).
    
    When used via the INFERNAL_PARALLEL subworkflow, receives individual chunks.
    Can also be called directly with full assemblies for single-node runs.
========================================================================================
*/

process INFERNAL {
    tag "${assembly_id}:${chunk_id}"
    label 'infernal'

    input:
    tuple val(assembly_id), val(chunk_id), path(masked_chunk)
    path rfam_database

    output:
    tuple val(assembly_id), val(chunk_id), path("${assembly_id}.${chunk_id}_rfam.tblout"), emit: infernal_tblout
    tuple val(assembly_id), val(chunk_id), path("${assembly_id}.${chunk_id}_rfam.cmscan"),  emit: infernal_cmscan

    script:
    """
    cmscan \
        --cpu ${task.cpus} \
        --cut_ga \
        --rfam \
        --nohmmonly \
        --tblout ${assembly_id}.${chunk_id}_rfam.tblout \
        --fmt 2 \
        --clanin ${rfam_database}/Rfam.clanin \
        ${rfam_database}/Rfam.cm \
        ${masked_chunk} > ${assembly_id}.${chunk_id}_rfam.cmscan
    """

    stub:
    """
    touch ${assembly_id}.${chunk_id}_rfam.tblout
    touch ${assembly_id}.${chunk_id}_rfam.cmscan
    """
}

/*
========================================================================================
    MERGE_INFERNAL MODULE
========================================================================================
    Merge Infernal results from multiple chunks into single files.

    tblout: Header lines start with '#'. Keep header from first file,
            concatenate data lines from all files.
    cmscan: Plain text alignment output — simple concatenation.
========================================================================================
*/

process MERGE_INFERNAL {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/ncRNA", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(tblout_files), path(cmscan_files)

    output:
    tuple val(assembly_id), path("${assembly_id}_rfam.tblout"), emit: infernal_tblout
    tuple val(assembly_id), path("${assembly_id}_rfam.cmscan"),  emit: infernal_cmscan

    script:
    """
    set -euo pipefail

    # ── Merge tblout files ────────────────────────────────────────────────────
    first_tblout=\$(ls -1 ${tblout_files} | sort -V | head -1)
    grep '^#' "\$first_tblout" > ${assembly_id}_rfam.tblout

    ls -1 ${tblout_files} | sort -V | while read -r f; do
        grep -v '^#' "\$f" >> ${assembly_id}_rfam.tblout || true
    done

    # ── Merge cmscan output ───────────────────────────────────────────────────
    > ${assembly_id}_rfam.cmscan
    ls -1 ${cmscan_files} | sort -V | while read -r f; do
        cat "\$f" >> ${assembly_id}_rfam.cmscan
    done
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