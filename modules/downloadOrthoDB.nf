/*
========================================================================================
    DOWNLOAD OrthoDB PROTEIN DATABASE
========================================================================================
    Downloads the specified OrthoDB clade partition for use with BRAKER.
    Pre-partitioned files are provided by the Greifswald team from OrthoDB v12.

    Available clades: Metazoa, Vertebrata, Viridiplantae, Arthropoda,
                      Fungi, Alveolata, Stramenopiles, Amoebozoa, Euglenozoa,
                      Eukaryota

    Stored permanently via storeDir — reused across pipeline runs.

    Inputs:
        - clade:  OrthoDB clade name (e.g. 'Vertebrata')

    Outputs:
        - orthodb_proteins:  ${clade}.fa

    Reference:
        Kuznetsov et al. (2023) OrthoDB v11. NAR 51(D1):D445-D451
        Pre-partitioned files: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/
========================================================================================
*/

process DOWNLOAD_ORTHODB {
    tag "${clade}"
    label 'download_db'
    storeDir "${params.orthodb_database}"

    input:
    val clade

    output:
    path("${clade}.fa"), emit: orthodb_proteins

    script:
    """
    set -euo pipefail

    wget --tries=10 --timeout=60 \\
        https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/${clade}.fa.gz

    gunzip ${clade}.fa.gz

    n_seqs=\$(grep -c '^>' ${clade}.fa)
    echo "Downloaded OrthoDB v12 ${clade}: \${n_seqs} protein sequences" >&2
    """

    stub:
    """
    touch ${clade}.fa
    """
}