/*
========================================================================================
    PREPARE PROTEIN DATABASE
========================================================================================
    Optionally download a reference proteome by NCBI accession and combine it
    with the OrthoDB protein database for BRAKER input.

    Handles four scenarios:
      1. OrthoDB + reference accession → combined database
      2. OrthoDB only                  → pass through OrthoDB
      3. Reference accession only      → download and use reference proteome
      4. Neither                       → produce NO_PROT_DB sentinel file

    The combined database is built in the work directory (not stored
    permanently) since the reference proteome is project-specific.

    Inputs:
        - orthodb_proteins:  OrthoDB clade FASTA (from DOWNLOAD_ORTHODB),
                             or NO_PROT_DB sentinel if no clade specified

    Outputs:
        - combined_proteins: Combined protein FASTA for BRAKER
                             (or NO_PROT_DB sentinel if no proteins at all)

    Tools: NCBI datasets CLI (conda: ncbi-datasets-cli)
========================================================================================
*/

process PREPARE_PROTEIN_DB {
    tag "${ref_accession ?: 'no_ref'}:${has_orthodb ? 'orthodb' : 'no_orthodb'}"
    label 'prepare_protein_db'

    input:
    path orthodb_proteins

    output:
    path("combined_proteins.fa"), emit: combined_proteins

    script:
    def ref_accession = params.ref_proteome_accession ?: ''
    // Check if orthodb input is a real file or the NO_PROT_DB sentinel
    has_orthodb = orthodb_proteins.name != 'NO_PROT_DB'

    if (ref_accession && has_orthodb)
        // Case 1: OrthoDB + reference proteome
        """
        set -euo pipefail

        # Download reference proteome
        datasets download genome accession ${ref_accession} \\
            --include protein \\
            --filename ref_proteome.zip

        unzip -o ref_proteome.zip -d ref_download

        ref_fasta=\$(find ref_download -name "protein.faa" -type f | head -1)
        if [ -z "\${ref_fasta}" ]; then
            echo "ERROR: No protein.faa found for accession ${ref_accession}" >&2
            exit 1
        fi

        n_ref=\$(grep -c '^>' "\${ref_fasta}")
        n_odb=\$(grep -c '^>' ${orthodb_proteins})
        echo "OrthoDB: \${n_odb}, Reference (${ref_accession}): \${n_ref}" >&2

        cat ${orthodb_proteins} "\${ref_fasta}" > combined_proteins.fa

        n_total=\$(grep -c '^>' combined_proteins.fa)
        echo "Combined: \${n_total} total proteins" >&2

        rm -rf ref_download ref_proteome.zip
        """

    else if (ref_accession && !has_orthodb)
        // Case 2: Reference proteome only (no OrthoDB)
        """
        set -euo pipefail

        datasets download genome accession ${ref_accession} \\
            --include protein \\
            --filename ref_proteome.zip

        unzip -o ref_proteome.zip -d ref_download

        ref_fasta=\$(find ref_download -name "protein.faa" -type f | head -1)
        if [ -z "\${ref_fasta}" ]; then
            echo "ERROR: No protein.faa found for accession ${ref_accession}" >&2
            exit 1
        fi

        cp "\${ref_fasta}" combined_proteins.fa

        n_ref=\$(grep -c '^>' combined_proteins.fa)
        echo "Reference proteome only (${ref_accession}): \${n_ref} proteins" >&2

        rm -rf ref_download ref_proteome.zip
        """

    else if (!ref_accession && has_orthodb)
        // Case 3: OrthoDB only
        """
        cp ${orthodb_proteins} combined_proteins.fa
        echo "Using OrthoDB proteins only (no reference proteome specified)" >&2
        """

    else
        // Case 4: No proteins at all
        """
        touch combined_proteins.fa
        echo "WARNING: No protein database specified. BRAKER will run without protein evidence." >&2
        """

    stub:
    """
    touch combined_proteins.fa
    """
}