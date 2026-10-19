/*
========================================================================================
    InterProScan MODULE
========================================================================================
    Run InterProScan on protein sequences
========================================================================================
*/

process INTERPROSCAN {
    tag "${assembly_id}:${chunk_id}"
    label 'interproscan'
    //publishDir "${params.outdir}/functional_annotation/${assembly_id}", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), val(chunk_id), path(protein_fasta)
    tuple val(database_source), path(interproscan_db)

    output:
    tuple val(assembly_id), val(chunk_id), path("${assembly_id}.${chunk_id}.interproscan.tsv"), emit: interproscan_tsv

    script:
    """
    set -euo pipefail

    test -d '${interproscan_db}/pfam'

    /opt/interproscan/interproscan.sh \\
        --input "${protein_fasta}" \\
        --output-file-base "${assembly_id}.${chunk_id}.interproscan" \\
        --cpu ${task.cpus} \\
        --formats TSV \\
        --disable-precalc \\
        --seqtype p \\
        --pathways \\
        --goterms \\
        --iprlookup \\
        --tempdir . \\
        --verbose
    """

    stub:
    """
    touch ${assembly_id}.${chunk_id}.interproscan.tsv
    """
}

/*
========================================================================================
    MERGE_INTERPROSCAN MODULE
========================================================================================
    Merge InterProScan results from multiple chunks into a single file
========================================================================================
*/

process MERGE_INTERPROSCAN {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/functional_annotation/${assembly_id}", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(tsvs)

    output:
    tuple val(assembly_id), path("${assembly_id}.interproscan.tsv"), emit: merged_tsv

    script:
    """
    set -euo pipefail
    cat ${tsvs.join(' ')} > ${assembly_id}.interproscan.tsv
    """
}

/*
========================================================================================
    DOWNLOAD InterProScan DATABASE MODULE
========================================================================================
    Downloads the specified InterProScan lineage database once, before parallel InterProScan runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_INTERPROSCAN_DB {
    tag "${interproscan_version}"
    label 'download_db'

    storeDir "${params.interproscan_database}"

    input:
    val interproscan_version

    output:
    path("interproscan-data-${interproscan_version}/data"), emit: db

    script:
    """
    set -euo pipefail

    url="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${interproscan_version}/alt/interproscan-data-${interproscan_version}.tar.gz"
    tarball="interproscan-data-${interproscan_version}.tar.gz"

    # download (quiet progress meter; retries; resume if partial exists)
    curl -fL -C - \\
      --retry 10 --retry-delay 10 --retry-all-errors \\
      --connect-timeout 30 --max-time 0 \\
      --silent --show-error \\
      -o "\$tarball" "\$url"

    # validate archive is readable
    tar -tzf "\$tarball" >/dev/null

    # extract
    tar -xzf "\$tarball"

    # find the extracted top-level dir WITHOUT using a pipe (no SIGPIPE possible)
    topdir=\$(find . -mindepth 1 -maxdepth 1 -type d -printf '%f\\n' -quit)
    test -n "\$topdir"

    # normalize folder name so Nextflow output path is deterministic
    if [ "\$topdir" != "interproscan-data-${interproscan_version}" ]; then
        mv "\$topdir" "interproscan-data-${interproscan_version}"
    fi

    # sanity check matches declared Nextflow output
    test -d "interproscan-data-${interproscan_version}/data"
    """
}
