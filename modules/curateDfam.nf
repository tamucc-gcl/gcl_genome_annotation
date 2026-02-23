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
    //publishDir "${params.outdir}/dfam", mode: params.publish_dir_mode

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

/*
========================================================================================
    DOWNLOAD DFam DATABASE MODULE
========================================================================================
*/

process DOWNLOAD_DFAM_DB {
  tag "${dfam_version}"
  label 'download_db'
  storeDir "${params.dfam_downloads}"

  input:
  val dfam_version

  output:
  path("dfam_db"), emit: dfam_db_dir

  script:
  """
  set -euo pipefail
  mkdir -p dfam_db
  cd dfam_db
  wget -r -np -nH --cut-dirs=4 -A "dfam${dfam_version}_full.*.h5.gz,dfam${dfam_version}_full.*.h5.gz.md5,README.txt" \
    https://www.dfam.org/releases/current/families/FamDB/

  gunzip dfam${dfam_version}_full.*.h5.gz
  """
}