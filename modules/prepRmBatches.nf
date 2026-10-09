// Preserve whole contigs so RepeatMasker never sees artificial tile boundaries.
process PREP_RM_BATCHES {
    tag "${assembly_id}"
    label 'prep_rm_batches'
    publishDir "${params.outdir}/repeats/batches", mode: 'copy', overwrite: true, pattern: '*.batches.json'
    input:
    tuple val(assembly_id), path(assembly_path)
    val n_chunks
    path helper
    output:
    tuple val(assembly_id), path('batches/batch_*.fa'), emit: batches
    tuple val(assembly_id), path('beds/batch_*.bed'), emit: beds
    tuple val(assembly_id), path("${assembly_id}_seq_lengths.tsv"), emit: seq_lengths
    tuple val(assembly_id), path("${assembly_id}.batches.json"), emit: accounting
    script:
    """
    python ${helper} batches --source '${assembly_path}' --assembly ${assembly_id} --count ${n_chunks} --outdir .
    """
}
