/*
========================================================================================
    AGAT MODULE
========================================================================================
    Get descriptive stats of structural annotations
========================================================================================
*/

process AGAT {
    tag "${assembly_id}"
    label 'agat'
    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), path(assembly_path), path(gff)
    
    output:
    tuple val(assembly_id), path("${assembly_id}_annotation_stats"), emit: annotation_stats
    tuple val(assembly_id), path("annotation_stats_distribution_plots/${assembly_id}"), emit: annotation_plots
    
    script:
    """
    grep -P "\\t(gene|mRNA|exon|CDS|intron|start_codon|stop_codon|tRNA|rRNA)\\t" ${gff} \
    > genes_only.gff3

    agat_sp_statistics.pl \
        --gff genes_only.gff3 \
        -f ${assembly_path} \
        -d \
        --yaml \
        --output ${assembly_id}_annotation_stats
    mkdir -p annotation_stats_distribution_plots
    mv ${assembly_id}_annotation_stats_distribution_plots annotation_stats_distribution_plots/${assembly_id}
    """
    
    stub:
    """
    touch ${assembly_id}_annotation_stats
    mkdir -p annotation_stats_distribution_plots/${assembly_id}
    """
}