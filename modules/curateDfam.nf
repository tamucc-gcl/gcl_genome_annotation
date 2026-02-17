/*
========================================================================================
    CURATE DFAM
========================================================================================
    Extract taxon-specific repeat library from Dfam database
    Runs ONCE per pipeline execution (not per genome)
    
    Inputs:
        - dfam_dir:   Path to Dfam database directory
        - taxa_name:  Taxonomic group name or NCBI taxid (e.g. 'Actinopterygii')
    
    Outputs:
        - dfam_lib: dfam_subset.lib (FASTA-format repeat library)
========================================================================================
*/

process CURATE_DFAM {
    label 'dfam'
    publishDir "${params.outdir}/dfam", mode: params.publish_dir_mode

    input:
    path dfam_dir
    val taxa_name

    output:
    path "dfam_subset.lib", emit: dfam_lib

    script:
    """
    famdb.py \
        -i ${dfam_dir} families \
        --format fasta_name \
        --include-class-in-name \
        --ancestors --descendants \
        ${taxa_name} > dfam_subset.lib
    """

    stub:
    """
    touch dfam_subset.lib
    """
}