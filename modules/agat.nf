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
    // Preserve the standalone UTR pilot's layout without a global selector for
    // its development-only alias. Main uses the final/selected QC scopes.
    publishDir "${params.outdir}/qc/${task.ext.qc_stage ?: (task.process == 'UTR_AGAT' ? 'utr' : 'baseline')}/${assembly_id}/agat", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), path(assembly_path), path(gff)
    
    output:
    tuple val(assembly_id), path("${assembly_id}_annotation_stats"), emit: annotation_stats
    tuple val(assembly_id), path("annotation_stats_distribution_plots/${assembly_id}"), emit: annotation_plots
    tuple val(assembly_id), path("${assembly_id}_annotation_stats.yaml"), emit: annotation_yaml
    
    script:
    """
    set -euo pipefail
    grep -P "\\t(gene|mRNA|exon|CDS|five_prime_UTR|three_prime_UTR|intron|start_codon|stop_codon|tRNA|rRNA)\\t" ${gff} \
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
    touch ${assembly_id}_annotation_stats.yaml
    mkdir -p annotation_stats_distribution_plots/${assembly_id}
    """
}
