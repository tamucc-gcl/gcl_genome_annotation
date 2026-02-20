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
    //publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), val(chunk_id), path(protein_fasta)
    path(interproscan_db)

    output:
    tuple val(assembly_id), val(chunk_id), path("${assembly_id}.${chunk_id}.interproscan.tsv"), emit: interproscan_tsv

    script:
    """
    set -euo pipefail

    # Strip '*' stop characters (InterProScan rejects them)
    awk '
      /^>/ {print; next}
      {gsub(/\\*/, "", \$0); print}
    ' "${protein_fasta}" > clean.nostop.fa

    /opt/interproscan/interproscan.sh \\
        --input clean.nostop.fa \\
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
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

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
    label 'interproscan_download'

    storeDir "${params.interproscan_downloads}"

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

    # extract to scratch
    tar -xzf "\$tarball"

    # normalize top-level folder name deterministically without SIGPIPE
    # (avoid `tar -t | head` under pipefail)
    topdir=\$(tar -tzf "\$tarball" | awk -F/ 'NF{print \$1; exit}')
    test -n "\$topdir"
    mv "\$topdir" "interproscan-data-${interproscan_version}"

    # sanity check matches declared Nextflow output
    test -d "interproscan-data-${interproscan_version}/data"
    """
}