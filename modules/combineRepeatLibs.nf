/*
========================================================================================
    COMBINE REPEAT LIBRARIES
========================================================================================
    Concatenate de novo curated repeats with Dfam taxon-specific library into a
    single combined library for RepeatMasker. Runs once per assembly rather than
    per batch.

    Inputs:
        - assembly_id:      Sample identifier
        - curated_repeats:  De novo curated repeat library
        - dfam_lib:         Dfam taxon-specific repeat library

    Outputs:
        - combined_lib:     Combined repeat library FASTA

    Tools: cat
========================================================================================
*/

process COMBINE_REPEAT_LIBS {
    tag "${assembly_id}"
    label 'combine_libs'
    publishDir "${params.outdir}/repeats/libraries", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(curated_repeats), path(dfam_lib)

    output:
    tuple val(assembly_id), path("${assembly_id}_combined_repeats.fa"), emit: combined_lib

    script:
    """
    set -euo pipefail
    cat ${curated_repeats} ${dfam_lib} > ${assembly_id}_combined_repeats.fa
    """

    stub:
    """
    touch ${assembly_id}_combined_repeats.fa
    """
}
