/*
========================================================================================
    DeepSIG MODULE
========================================================================================
    Run DeepSIG on protein sequences
========================================================================================
*/

process DEEPSIG {
    tag "${assembly_id}:${chunk_id}"
    label 'deepsig'
    //publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), val(chunk_id), path(protein_fasta)

    output:
    tuple val(assembly_id), val(chunk_id),
        path("${assembly_id}_chunk_${chunk_id}_deepsig_annotation.gff3"),
        emit: deepsig_annotations

    script:
    """
    set -euo pipefail

    deepsig \
        --fasta ${protein_fasta} \
        --outf ${assembly_id}_chunk_${chunk_id}_deepsig_annotation.gff3 \
        --organism euk \
        --outfmt gff3 \
        --threads ${task.cpus}
    """

    stub:
    """
    touch ${assembly_id}_chunk_${chunk_id}_deepsig_annotation.gff3
    """
}

/*
========================================================================================
    MERGE_DEEPSIG MODULE
========================================================================================
    Merge DeepSIG results from multiple chunks into a single file
========================================================================================
*/

process MERGE_DEEPSIG {
    tag "${assembly_id}"
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(gff_files)

    output:
    tuple val(assembly_id), path("${assembly_id}.deepsig.gff3"), emit: merged_gff3

    script:
    """
    set -euo pipefail
    
    ls -1 ${gff_files} | sort -V > gff.list

    first=\$(head -n 1 gff.list)

    # Write header lines from first chunk
    grep '^#' "\$first" > ${assembly_id}.deepsig.gff3

    # Append data lines (skip comment/header lines) from ALL chunks including first
    while read -r f; do
        grep -v '^#' "\$f" >> ${assembly_id}.deepsig.gff3 || true
    done < gff.list
    """
}
