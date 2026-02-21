/*
========================================================================================
    eggNOG MODULE
========================================================================================
    Run eggNOG on protein sequences
========================================================================================
*/

process EGGNOG {
  tag "${assembly_id}:${chunk_id}"
  label 'eggnog'
  //publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

  input:
  tuple val(assembly_id), val(chunk_id), path(protein_fasta)
  path eggnog_db_dir

  output:
  tuple val(assembly_id), val(chunk_id),
    path("${assembly_id}_chunk_${chunk_id}_eggnog.emapper.annotations"),
    emit: eggnog_annotations

  script:
  """
  set -euo pipefail

  emapper.py \
    -i "${protein_fasta}" \\
    --itype proteins \\
    -o "${assembly_id}_chunk_${chunk_id}_eggnog" \\
    --cpu ${task.cpus} \\
    --data_dir "${eggnog_db_dir}" \\
    -m diamond \\
    --evalue "${params.eggnog_evalue}" \\
    --sensmode "${params.eggnog_sensemode}" \\
    --dmnd_iterate yes \\
    --dmnd_algo auto \\
    --pfam_realign none \\
    --tax_scope "${params.eggnog_taxid}" \\
    --target_orthologs "${params.eggnog_target_orthologs}" \\
    --go_evidence "${params.eggnog_goEvidence}"
  """

  stub:
  """
  touch ${assembly_id}_chunk_${chunk_id}_eggnog.emapper.annotations
  """
}

/*
========================================================================================
    MERGE_EGGNOG MODULE
========================================================================================
    Merge eggNOG results from multiple chunks into a single file
========================================================================================
*/

process MERGE_EGGNOG {
    tag "${assembly_id}"
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(annot_files)

    output:
    path("${assembly_id}.eggnog.annotations"), emit: merged

    script:
    """
    set -euo pipefail

    # Ensure deterministic order (part_001, part_002, ...)
    ls -1 ${annot_files} | sort -V > annot.list

    first=\$(head -n 1 annot.list)

    # Write header + first chunk's data
    cat "\$first" > ${assembly_id}.eggnog.annotations

    # Append ONLY data rows (skip all lines beginning with '#') from the remaining chunks
    tail -n +2 annot.list | while read -r f; do
        awk 'BEGIN{FS="\\t"} !/^#/' "\$f" >> ${assembly_id}.eggnog.annotations
    done

    # Optional sanity: remove any accidental blank lines
    sed -i '/^\\s*\$/d' ${assembly_id}.eggnog.annotations
    """
}

/*
========================================================================================
    DOWNLOAD eggNOG DATABASE MODULE
========================================================================================
    Downloads the specified eggNOG database once, before parallel eggNOG runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_EGGNOG_DB {
  tag "${eggnog_version}"
  label 'download_db'
  storeDir "${params.eggnog_downloads}"

  input:
  val eggnog_version

  output:
  path("eggnog_db"), emit: eggnog_db_dir

  script:
  """
  set -euo pipefail

  # Download directly into storeDir (i.e. this directory)
  wget --tries=10 --timeout=30 -c -O eggnog.db.gz            http://eggnog5.embl.de/download/emapperdb-${eggnog_version}/eggnog.db.gz
  wget --tries=10 --timeout=30 -c -O eggnog.taxa.tar.gz      http://eggnog5.embl.de/download/emapperdb-${eggnog_version}/eggnog.taxa.tar.gz
  wget --tries=10 --timeout=30 -c -O eggnog_proteins.dmnd.gz http://eggnog5.embl.de/download/emapperdb-${eggnog_version}/eggnog_proteins.dmnd.gz

  gunzip -f eggnog.db.gz
  tar -xzf eggnog.taxa.tar.gz && rm -f eggnog.taxa.tar.gz
  gunzip -f eggnog_proteins.dmnd.gz

  # Emit a named directory path without creating a subdirectory:
  # eggnog_db is a symlink pointing to "this dir"
  ln -sfn . eggnog_db
  """
}