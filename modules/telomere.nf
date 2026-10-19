// Reuse tidk search/plot; conversion follows native end-coordinate semantics.
process TELOMERE_ANNOTATION {
    tag "${assembly_id}"
    label 'telomere'
    publishDir "${params.outdir}/telomeres", mode: params.publish_dir_mode
    input:
    tuple val(assembly_id), path(assembly_fasta), val(motif)
    path helpers
    output:
    tuple val(assembly_id), path("${assembly_id}.telomeres.gff3"), emit: telomere_gff3
    tuple val(assembly_id), path("${assembly_id}.telomeres.tsv"), emit: telomere_tsv
    tuple val(assembly_id), path("${assembly_id}.telomeres.svg"), emit: telomere_plot
    script:
    """
    set -euo pipefail
    tidk search --string '${motif}' --window ${params.telomere_window} --output '${assembly_id}' --dir tidk_out '${assembly_fasta}'
    cp 'tidk_out/${assembly_id}_telomeric_repeat_windows.tsv' '${assembly_id}.telomeres.tsv'
    tidk plot --tsv '${assembly_id}.telomeres.tsv' --output '${assembly_id}.telomeres'
    python genome_tracks.py telomeres --tsv '${assembly_id}.telomeres.tsv' --genome '${assembly_fasta}' --aid '${assembly_id}' --motif '${motif}' --window ${params.telomere_window} --minimum ${params.telomere_min_count} --terminal ${params.telomere_end_window} --output '${assembly_id}.telomeres.gff3'
    """
    stub:
    """
    touch ${assembly_id}.telomeres.gff3 ${assembly_id}.telomeres.tsv ${assembly_id}.telomeres.svg
    """
}
