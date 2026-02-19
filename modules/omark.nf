/*
========================================================================================
    OMArk MODULE
========================================================================================
    Assess genome completeness using OMArk
========================================================================================
*/

process OMARK {
    tag "${assembly_id}"
    label 'omark'
    publishDir "${params.outdir}/qc/omark", mode: params.publish_dir_mode
    
    input:
    tuple val(assembly_id), path(protein_fasta)
    path(omamer_db)
    path(taxa_db)
    
    output:
    tuple val(assembly_id), path("${assembly_id}.png"), emit: omark_plot
    tuple val(assembly_id), path("${assembly_id}_detailed_summary.txt"), emit: omark_detailed
    tuple val(assembly_id), path("${assembly_id}.sum"), emit: omark_summary
    
    script:
    """
    omamer search \
        --nthreads ${task.cpus} \\
        --db ${omamer_db} \\
        --query ${protein_fasta} \\
        --out ${assembly_id}.omamer

    omark \
        --file ${assembly_id}.omamer \\
        --database ${omamer_db} \\
        --ete_ncbi_db ${taxa_db} \\
        --taxid ${params.ncbi_speciesID} \\
        --og_fasta ${protein_fasta}
    
    mv omark_output/${assembly_id}.png ${assembly_id}.png
    mv omark_output/${assembly_id}_detailed_summary.txt ${assembly_id}_detailed_summary.txt
    mv omark_output/${assembly_id}.sum ${assembly_id}.sum
    """
    
    stub:
    """
    touch ${assembly_id}.png
    touch ${assembly_id}_detailed_summary.txt
    touch ${assembly_id}.sum
    """
}