/*
========================================================================================
    CURATE REPEATS
========================================================================================
    Filter, deduplicate, cluster, and classify de novo repeat families
    
    Inputs:
        - assembly_id:      Sample identifier
        - families_fa:      ${assembly_id}-families.fa from RepeatModeler
        - cluster_pct:      CD-HIT-EST clustering threshold (0-1); null/0 = skip clustering
        - min_repeat_length: Minimum repeat length to retain (default: 100)
    
    Outputs:
        - curated_repeats: ${assembly_id}_curatedRepeats.fasta
========================================================================================
*/

process CURATE_REPEATS {
    tag "${assembly_id}"
    label 'repeat'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(families_fa)

    output:
    tuple val(assembly_id), path("${assembly_id}_curatedRepeats.fasta"), emit: curated_repeats

    script:
    def min_len = params.min_repeat_length ?: 100
    def cluster_pct = params.cluster_pct ?: 0
    def do_cluster = cluster_pct > 0 && cluster_pct <= 1
    """
    # 1. Remove short families
    seqkit seq -m ${min_len} -g ${families_fa} > rm_filtered_len.fa

    # 2. Remove simple & low complexity repeats
    grep -v -E "Simple_repeat|Low_complexity|Satellite" rm_filtered_len.fa > rm_filtered_nosimple.fa || true

    # 3. Remove exact duplicates
    seqkit rmdup -s rm_filtered_nosimple.fa > rm_filtered_nodup.fa

    # 4. Cluster near-identical families (optional)
    if [ "${do_cluster}" = "true" ]; then
        cd-hit-est \
            -i rm_filtered_nodup.fa \
            -o rm_clustered.fa \
            -c ${cluster_pct} \
            -T ${task.cpus} \
            -M 0
    else
        cp rm_filtered_nodup.fa rm_clustered.fa
    fi

    # 5. Reclassify with RepeatClassifier
    RepeatClassifier -consensi rm_clustered.fa

    cp rm_clustered.fa.classified ${assembly_id}_curatedRepeats.fasta
    """

    stub:
    """
    touch ${assembly_id}_curatedRepeats.fasta
    """
}