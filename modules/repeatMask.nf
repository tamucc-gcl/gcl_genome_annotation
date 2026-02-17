/*
========================================================================================
    REPEAT MASKER
========================================================================================
    Mask repeats in genome assembly using combined de novo + Dfam libraries
    
    Inputs:
        - assembly_id:      Sample identifier
        - assembly_path:    Path to genome assembly FASTA
        - curated_repeats:  De novo curated repeat library
        - dfam_lib:         Dfam taxon-specific repeat library
        - job_speed:        RepeatMasker speed flag (null, 's', 'q', 'qq')
    
    Outputs:
        - masked_fasta: ${assembly_id}.masked.fasta (soft-masked genome)
        - repeat_gff:   ${assembly_id}.gff          (repeat annotation)
        - repeat_tbl:   ${assembly_id}.tbl           (repeat summary table)
========================================================================================
*/

process REPEAT_MASK {
    tag "${assembly_id}"
    label 'repeat'
    publishDir "${params.outdir}/repeats", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(assembly_path), path(curated_repeats)
    path dfam_lib

    output:
    tuple val(assembly_id), path("${assembly_id}.masked.fasta"), emit: masked_fasta
    tuple val(assembly_id), path("${assembly_id}.gff"),          emit: repeat_gff
    tuple val(assembly_id), path("${assembly_id}.tbl"),          emit: repeat_tbl

    script:
    def speed_flag = params.job_speed ? "-${params.job_speed}" : "-s"
    def rm_threads = Math.max(1, (task.cpus / 4) as int)
    """
    # 1. Combine de novo + Dfam repeat libraries
    cat ${curated_repeats} ${dfam_lib} > combined_repeats.fa

    # 2. Run RepeatMasker
    RepeatMasker \
        -parallel ${rm_threads} \
        -engine ncbi \
        ${speed_flag} \
        -cutoff 225 \
        -lib combined_repeats.fa \
        -xsmall -small \
        -gff \
        ${assembly_path}

    # 3. Collect outputs
    cp ${assembly_path}.out.gff ${assembly_id}.gff
    cp ${assembly_path}.masked  ${assembly_id}.masked.fasta
    cp ${assembly_path}.tbl     ${assembly_id}.tbl
    """

    stub:
    """
    touch ${assembly_id}.masked.fasta
    touch ${assembly_id}.gff
    touch ${assembly_id}.tbl
    """
}