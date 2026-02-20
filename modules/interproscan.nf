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
    publishDir "${params.outdir}/qc/interproscan", mode: params.publish_dir_mode

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

process MERGE_INTERPROSCAN {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/qc/interproscan", mode: params.publish_dir_mode

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