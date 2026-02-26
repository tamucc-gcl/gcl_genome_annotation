/*
========================================================================================
    FILTER ASSEMBLY
========================================================================================
    Filter genome assembly contigs by minimum length and report basic assembly
    statistics of the filtered output.

    Inputs:
        - assembly_id:       Sample identifier
        - assembly_path:     Path to genome assembly FASTA
        - min_contig_length: Minimum contig length to retain (val, from params)

    Outputs:
        - filtered_assembly: ${assembly_id}.filtered.fasta  (contigs >= min length)
        - assembly_stats:    ${assembly_id}.assembly_stats.tsv (basic stats)

    Tools: seqkit (conda: bioconda::seqkit), assembly-stats (conda: bioconda::assembly-stats)
========================================================================================
*/

process FILTER_ASSEMBLY {
    tag "${assembly_id}"
    label 'filter_assembly'
    publishDir "${params.outdir}/filtered_assemblies", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(assembly_path)

    output:
    tuple val(assembly_id), path("${assembly_id}.filtered.fasta"), emit: filtered_assembly
    tuple val(assembly_id), path("${assembly_id}.assembly_stats.tsv"), emit: assembly_stats

    script:
    """
    set -euo pipefail

    # ── 1. Filter contigs by minimum length ───────────────────────────────────
    seqkit seq -m ${params.min_contig_length} -g "${assembly_path}" > ${assembly_id}.filtered.fasta

    # ── 2. Count removed contigs for logging ──────────────────────────────────
    total_before=\$(grep -c '^>' "${assembly_path}")
    total_after=\$(grep -c '^>' ${assembly_id}.filtered.fasta)
    removed=\$((total_before - total_after))
    echo "Filtered ${assembly_id}: kept \${total_after}/\${total_before} contigs (removed \${removed} < ${params.min_contig_length} bp)" >&2

    # ── 3. Generate assembly stats ────────────────────────────────────────────
    assembly-stats -t ${assembly_id}.filtered.fasta > ${assembly_id}.assembly_stats.tsv
    """

    stub:
    """
    touch ${assembly_id}.filtered.fasta
    touch ${assembly_id}.assembly_stats.tsv
    """
}