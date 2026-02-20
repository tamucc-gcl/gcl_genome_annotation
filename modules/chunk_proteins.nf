/*
========================================================================================
    Chunk Proteins MODULE
========================================================================================
    Split protein sequences into smaller chunks for parallel processing in InterProScan, Diamond, eggNOG
========================================================================================
*/

process CHUNK_PROTEINS {
    tag "${assembly_id}"
    label 'chunk_proteins'
    //publishDir "${params.outdir}/chunk_proteins", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(protein_fasta)

    output:
    tuple val(assembly_id), path("${assembly_id}_chunks"), emit: chunk_dir
    tuple val(assembly_id), path("${assembly_id}_chunks/*"), emit: chunks
    tuple val(assembly_id), path("${assembly_id}_chunks.manifest.txt"), emit: manifest

    script:
    """
    mkdir -p ${assembly_id}_chunks
    seqkit split2 -j ${task.cpus} -p ${params.n_chunks} -O ${assembly_id}_chunks "${protein_fasta}"

    find ${assembly_id}_chunks -maxdepth 1 -type f \\( \
        -name "*.fa" -o -name "*.faa" -o -name "*.fasta" -o \
        -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.faa.gz" \
      \\) | sort > ${assembly_id}_chunks.manifest.txt
    """

    stub:
    """
    mkdir -p ${assembly_id}_chunks
    cat > ${assembly_id}_chunks/${assembly_id}.chunk_0001.faa << 'EOF'
    >fake_1
    MADEUP
    EOF
    echo "${assembly_id}_chunks/${assembly_id}.chunk_0001.faa" > ${assembly_id}_chunks.manifest.txt
    """
}