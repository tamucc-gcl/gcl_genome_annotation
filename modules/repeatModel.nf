/*
========================================================================================
    REPEAT MODELER
========================================================================================
    De novo repeat family identification using RepeatModeler2
    
    Inputs:
        - assembly_id:   Sample identifier
        - assembly_path: Path to genome assembly FASTA
    
    Outputs:
        - families_fa:  ${assembly_id}-families.fa  (consensus sequences per family)
        - families_stk: ${assembly_id}-families.stk (seed alignments per family)
        - log_file:     ${assembly_id}-rmod.log     (execution log)
========================================================================================
*/

process REPEAT_MODEL {
    tag "${assembly_id}"
    label 'repeat'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    container "${params.repeatmodeler_container}"

    input:
    tuple val(assembly_id), path(assembly_path)

    output:
    tuple val(assembly_id), path("${assembly_id}-families.fa"),  emit: families_fa
    tuple val(assembly_id), path("${assembly_id}-families.stk"), emit: families_stk
    tuple val(assembly_id), path("${assembly_id}-rmod.log"),     emit: log_file

    script:
    """
    export BLAST_USAGE_REPORT=false

    # Build BLAST database for RepeatModeler
    BuildDatabase \
        -name ${assembly_id} \
        ${assembly_path}

    # Run RepeatModeler with LTR structural search
    RepeatModeler \
        -threads ${task.cpus} \
        -database ${assembly_id} \
        -LTRStruct \
        -numAddlRounds ${params.repModel_rounds} 2>&1 | tee repeatmodeler.log

    # Collect outputs (RepeatModeler places them in RM_ directories)
    # cp RM_*/families-rmod.fa  ${assembly_id}-families.fa
    # cp RM_*/families-rmod.stk ${assembly_id}-families.stk
    # cp repeatmodeler.log       ${assembly_id}-rmod.log
    """

    stub:
    """
    touch ${assembly_id}-families.fa
    touch ${assembly_id}-families.stk
    touch ${assembly_id}-rmod.log
    """
}