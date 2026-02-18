/*
========================================================================================
    CURATE REPEATS
========================================================================================
    Filter, deduplicate, cluster, and classify de novo repeat families
    
    Inputs:
        - assembly_id:      Sample identifier
        - families_fa:      ${assembly_id}-families.fa from RepeatModeler
    
    Outputs:
        - curated_repeats: ${assembly_id}_curatedRepeats.fasta
========================================================================================
*/

process CLASSIFY_REPEATS {
    tag "${assembly_id}"
    label 'repeat_classify'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(families_fa)

    output:
    tuple val(assembly_id), path("${assembly_id}_curatedRepeats.fasta"), emit: curated_repeats

    script:
    """
    # Reclassify with RepeatClassifier
    RepeatClassifier -consensi ${families_fa}

    cp ${families_fa}.classified ${assembly_id}_curatedRepeats.fasta
    """

    stub:
    """
    touch ${assembly_id}_curatedRepeats.fasta
    """
}