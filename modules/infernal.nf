/*
========================================================================================
    Infernal
========================================================================================
    ncRNA prediction using Infernal with the Rfam database. This step identifies non-coding RNA genes in the genome assembly, including tRNAs, rRNAs, and other small RNAs.
    
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
    tRNAscan-SE \
        --thread ${task.cpus} \
        --eukaryotic \
        --isotype \
        --detail \
        --output ${assembly_id}.trnascan.tsv \
        --gff ${assembly_id}.trnascan.gff3 \
        --fasta ${assembly_id}.trnascan.fa \
        --bed ${assembly_id}.trnascan.bed \
        ${masked_assembly}
    """

    stub:
    """
    touch ${assembly_id}.trnascan.gff3
    touch ${assembly_id}.trnascan.fa
    touch ${assembly_id}.trnascan.bed
    touch ${assembly_id}.trnascan.tsv
    """
}