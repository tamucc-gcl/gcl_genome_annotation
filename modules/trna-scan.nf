/*
========================================================================================
    tRNA-Scan
========================================================================================
    tRNA prediction using tRNA-Scan-SE
    
    Inputs:
        - assembly_id:      Sample identifier
        - masked_assembly:  Soft-masked genome assembly FASTA
    
    Outputs:
        - trnascan_gff3:     ${assembly_id}.trnascan.gff3           (tRNA models)
        - trnascan_fa:       ${assembly_id}.trnascan.fa              (tRNA sequences)
        - trnascan_bed:      ${assembly_id}.trnascan.bed             (tRNA BED format)
        - trnascan_tsv:      ${assembly_id}.trnascan.tsv             (tRNA TSV format)
========================================================================================
*/

process tRNASCAN {
    tag "${assembly_id}"
    label 'trna_scan'
    publishDir "${params.outdir}/ncRNA", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(masked_assembly)

    output:
    tuple val(assembly_id), path("${assembly_id}.trnascan.gff3"),            emit: trnascan_gff3
    tuple val(assembly_id), path("${assembly_id}.trnascan.fa"),              emit: trnascan_fa
    tuple val(assembly_id), path("${assembly_id}.trnascan.bed"), emit: trnascan_bed
    tuple val(assembly_id), path("${assembly_id}.trnascan.tsv"), emit: trnascan_tsv

    script:
    """
    set -euo pipefail
    tRNAscan-SE \
        --thread ${task.cpus} \
        -E \
        --detail \
        --output ${assembly_id}.trnascan.tsv \
        --gff ${assembly_id}.trnascan.gff3 \
        --fasta ${assembly_id}.trnascan.fa \
        --bed ${assembly_id}.trnascan.bed \
        ${masked_assembly}
    # Native successful no-hit runs may omit these optional files.
    touch ${assembly_id}.trnascan.gff3 ${assembly_id}.trnascan.fa ${assembly_id}.trnascan.bed ${assembly_id}.trnascan.tsv
    """

    stub:
    """
    touch ${assembly_id}.trnascan.gff3
    touch ${assembly_id}.trnascan.fa
    touch ${assembly_id}.trnascan.bed
    touch ${assembly_id}.trnascan.tsv
    """
}