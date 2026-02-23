/*
========================================================================================
    Diamond MODULE
========================================================================================
    Run Diamond on protein sequences
========================================================================================
*/

process DIAMOND {
  tag "${assembly_id}:${chunk_id}"
  label 'diamond'
  //publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

  input:
  tuple val(assembly_id), val(chunk_id), path(protein_fasta)
  path diamond_db

  output:
  tuple val(assembly_id), val(chunk_id),
    path("${assembly_id}.${chunk_id}.diamond.tsv"),
    emit: diamond_annotations

  script:
  """
  set -euo pipefail

  diamond blastp \
    --threads ${task.cpus} \
    --db ${diamond_db} \
    --query ${protein_fasta} \
    --out ${assembly_id}.${chunk_id}.diamond.tsv \
    --outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --${params.diamond_sensitivity} \
    --max-target-seqs ${params.diamond_maxHits} \
    --evalue ${params.diamond_evalue} \
    --query-cover ${params.diamond_qcov}

  """

  stub:
  """
  touch ${assembly_id}.${chunk_id}.diamond.tsv
  """
}

/*
========================================================================================
    MERGE_DIAMOND MODULE
========================================================================================
    Merge diamond results from multiple chunks into a single file
========================================================================================
*/

process MERGE_DIAMOND {
    tag "${assembly_id}"
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(annot_files)

    output:
    path("${assembly_id}.diamond.tsv"), emit: merged

    script:
    """
    set -euo pipefail
    # Ensure deterministic order (part_001, part_002, ...)
    ls -1 ${annot_files} | sort -V > annot.list

    first=\$(head -n 1 annot.list)

    # Write header + first chunk's data
    cat "\$first" > ${assembly_id}.diamond.tsv

    # Append ONLY data rows (skip all lines beginning with '#') from the remaining chunks
    tail -n +2 annot.list | while read -r f; do
        awk 'BEGIN{FS="\\t"} !/^#/' "\$f" >> ${assembly_id}.diamond.tsv
    done

    # Optional sanity: remove any accidental blank lines
    sed -i '/^\\s*\$/d' ${assembly_id}.diamond.tsv
    """
}

/*
========================================================================================
    DOWNLOAD Diamond DATABASE MODULE
========================================================================================
    Downloads the specified Diamond database once, before parallel diamond runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_DIAMOND_DB {
  tag "${diamond_version}"
  label 'download_diamond_db'
  storeDir "${params.diamond_database}"

  input:
  val diamond_version

  output:
  path("${diamond_version}.dmnd"), emit: diamond_db

  script:
  """
  set -euo pipefail
  wget --tries=10 --timeout=30 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

  diamond makedb \\
    --threads ${task.cpus} \\
    --db ${diamond_version} \\
    --in uniprot_sprot.fasta.gz
  """
}