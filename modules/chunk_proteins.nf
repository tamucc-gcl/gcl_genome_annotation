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

  //publishDir "${params.outdir}/protein_chunks", mode: params.publish_dir_mode

  input:
  tuple val(assembly_id), path(protein_fasta)

  output:
  tuple val(assembly_id), path("${assembly_id}_chunks/*"), emit: chunks

  script:
  """
  mkdir -p ${assembly_id}_chunks
  seqkit split2 -j ${task.cpus} -p ${params.n_chunks} -O ${assembly_id}_chunks "${protein_fasta}"
  """
}