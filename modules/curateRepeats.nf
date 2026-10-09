/*
========================================================================================
    CURATE REPEATS
========================================================================================
    Filter, deduplicate and optionally cluster de novo repeat families
    
    Inputs:
        - assembly_id:      Sample identifier
        - families_fa:      ${assembly_id}-families.fa from RepeatModeler
        - cluster_pct:      CD-HIT-EST clustering threshold (0-1); null/0 = skip clustering
        - min_repeat_length: Minimum repeat length to retain (default: 100)
    
    Outputs:
        - deduplicated and curated repeats: ${assembly_id}_rm_clustered.fasta
========================================================================================
*/

process CURATE_REPEATS {
    tag "${assembly_id}"
    label 'curate_repeat'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(families_fa)
    path helper

    output:
    tuple val(assembly_id), path("${assembly_id}_rm_clustered.fasta"), emit: clustered_repeats
    tuple val(assembly_id), path("${assembly_id}_repeat_curation.tsv"), emit: decisions

    script:
    def min_len = params.min_repeat_length == null ? 100 : params.min_repeat_length
    def cluster_pct = params.cluster_pct ?: 0
    def do_cluster = cluster_pct > 0 && cluster_pct <= 1
    """
    set -euo pipefail
    # Filter complete FASTA records; never remove a header while retaining its sequence.
    python ${helper} curate --source '${families_fa}' --minimum ${min_len} \
        --output rm_filtered_nodup.fa --report ${assembly_id}_repeat_curation.tsv

    # 4. Cluster near-identical families (optional)
    if [ "${do_cluster}" = "true" ] && [ -s rm_filtered_nodup.fa ]; then
        cd-hit-est \
            -i rm_filtered_nodup.fa \
            -o ${assembly_id}_rm_clustered.fasta \
            -c ${cluster_pct} \
            -T ${task.cpus} \
            -M 0
    else
        cp rm_filtered_nodup.fa ${assembly_id}_rm_clustered.fasta
    fi

    """

    stub:
    """
    touch ${assembly_id}_rm_clustered.fasta
    touch ${assembly_id}_repeat_curation.tsv
    """
}
