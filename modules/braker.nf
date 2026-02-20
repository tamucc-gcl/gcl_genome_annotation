/*
========================================================================================
    BRAKER3
========================================================================================
    Gene prediction using BRAKER3 with protein evidence and AUGUSTUS
    
    Inputs:
        - assembly_id:      Sample identifier
        - masked_assembly:  Soft-masked genome assembly FASTA
        - protein_database: Protein sequences for evidence (e.g. OrthoDB Vertebrata.fa)
        - species_name:     Species name for AUGUSTUS (e.g. 'Spratelloides_delicatulus')
        - augustus_config:   Path to AUGUSTUS config directory
    
    Outputs:
        - braker_gff3:     ${assembly_id}_braker.gff3           (gene models)
        - braker_proteins: ${assembly_id}_braker.aa              (predicted proteins)
        - braker_cds:      ${assembly_id}_braker.codingseq.fasta (coding sequences)
========================================================================================
*/

process BRAKER {
    tag "${assembly_id}"
    label 'braker'
    publishDir "${params.outdir}/structural_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(masked_assembly)
    path protein_database
    //path augustus_config

    output:
    tuple val(assembly_id), path("${assembly_id}.gff3"),            emit: braker_gff3
    tuple val(assembly_id), path("${assembly_id}_proteins.fa"),              emit: braker_proteins
    tuple val(assembly_id), path("${assembly_id}_cds.fa"), emit: braker_cds

    script:
    def prot_arg       = protein_database.name != 'NO_PROT_DB' ? "--prot_seq=${protein_database}" : ''
    def ab_initio_arg  = params.braker_ab_initio ? '--AUGUSTUS_ab_initio' : ''
    def rounds         = params.braker_rounds ?: 5
    //def species        = params.braker_species ?: assembly_id
    """
    # Copy AUGUSTUS config from container to writable scratch location
    mkdir -p augustus_config
    cp -a /opt/Augustus/config/. augustus_config/
    chmod -R u+rwX augustus_config
    ls -a
    ls -a augustus_config

    braker.pl \
        --threads=${task.cpus} \
        --genome=${masked_assembly} \
        ${prot_arg} \
        --species=${assembly_id} \
        --gff3 \
        --AUGUSTUS_CONFIG_PATH=\${PWD}/augustus_config \
        --rounds=${rounds} \
        ${ab_initio_arg} \
        --verbosity=3

    mv braker/braker.gff3     ${assembly_id}.gff3
    mv braker/braker.aa       ${assembly_id}_proteins.fa
    mv braker/braker.codingseq ${assembly_id}_cds.fa
    """

    stub:
    """
    touch ${assembly_id}.gff3
    touch ${assembly_id}_proteins.fa
    touch ${assembly_id}_cds.fa
    """
}