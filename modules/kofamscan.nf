/*
========================================================================================
    KoFamScan MODULE
========================================================================================
    Run KoFamScan on protein sequences
========================================================================================
*/
process KOFAMSCAN {
    tag "${assembly_id}:${chunk_id}"
    label 'kofamscan'
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), val(chunk_id), path(protein_fasta)
    path kofamscan_db_dir

    output:
    tuple val(assembly_id), val(chunk_id),
        path("${assembly_id}_chunk_${chunk_id}_kofamscan.tsv"),
        emit: kofamscan_annotations

    script:
    """
    set -euo pipefail

    exec_annotation \
        --cpu ${task.cpus} \
        --format detail-tsv \
        --profile ${kofamscan_db_dir}/profiles \
        --ko-list ${kofamscan_db_dir}/ko_list \
        -o ${assembly_id}_chunk_${chunk_id}_kofamscan.tsv \
        ${protein_fasta}
    """

    stub:
    """
    touch ${assembly_id}_chunk_${chunk_id}_kofamscan.tsv
    """
}

/*
========================================================================================
    MERGE_KOFAMSCAN MODULE
========================================================================================
    Merge KoFamScan results from multiple chunks into a single file
========================================================================================
*/

process MERGE_KOFAMSCAN {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(annot_files)

    output:
    tuple val(assembly_id), path("${assembly_id}.kofamscan.tsv"), emit: merged

    script:
    """
    set -euo pipefail

    ls -1 ${annot_files} | sort -V > annot.list

    first=\$(head -n 1 annot.list)

    # Write header from first chunk
    head -n 1 "\$first" > ${assembly_id}.kofamscan.tsv

    # Append data rows from all chunks
    while read -r f; do
        tail -n +2 "\$f" >> ${assembly_id}.kofamscan.tsv
    done < annot.list

    sed -i '/^\\s*\$/d' ${assembly_id}.kofamscan.tsv
    """
}

/*
========================================================================================
    DOWNLOAD KoFamScan DATABASE MODULE
========================================================================================
    Downloads the specified KoFamScan database once, before parallel KoFamScan runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_KOFAMSCAN_DB {
    tag "${kofamscan_version}"
    label 'download_db'
    storeDir "${params.kofamscan_database}"

    input:
    val kofamscan_version

    output:
    path("kofamscan_db"), emit: kofamscan_db_dir

    script:
    """
    set -euo pipefail

    mkdir -p kofamscan_db
    cd kofamscan_db
    wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz

    gunzip ko_list.gz
    tar -xzf profiles.tar.gz
    rm -f profiles.tar.gz
    """
}