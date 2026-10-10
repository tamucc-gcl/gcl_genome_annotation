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
    publishDir "${params.outdir}/repeats/denovo/${assembly_id}", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(assembly_path), path(dfam_source)
    val seed

    output:
    tuple val(assembly_id), path("${assembly_id}-families.fa"),  emit: families_fa
    tuple val(assembly_id), path("${assembly_id}-families.stk"), emit: families_stk
    tuple val(assembly_id), path("${assembly_id}-rmod.log"),     emit: log_file
    tuple val(assembly_id), path('RM_*', optional: true), emit: training
    tuple val(assembly_id), path('repeatmodeler.log'), emit: execution_log

    script:
    """
    set -euo pipefail
    export BLAST_USAGE_REPORT=false

    # Build BLAST database for RepeatModeler
    BuildDatabase \
        -name ${assembly_id} \
        ${assembly_path}

    # Run RepeatModeler with LTR structural search
    RepeatModeler \
        -threads ${task.cpus} \
        -database ${assembly_id} \
        -srand ${seed} \
        -LTRStruct \
        -numAddlRounds ${params.repModel_rounds} 2>&1 | tee repeatmodeler.log

    # RepeatModeler writes the declared families/Stockholm/log outputs directly.
    # Retain its RM_* training directories through Nextflow's native stage-out.
    """

    stub:
    """
    touch ${assembly_id}-families.fa
    touch ${assembly_id}-families.stk
    touch ${assembly_id}-rmod.log
    touch repeatmodeler.log
    """
}
