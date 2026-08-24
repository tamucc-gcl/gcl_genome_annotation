/*
========================================================================================
    PREP REPEATMASKER BATCHES
========================================================================================
    Split a genome assembly into N equal-sized chunks for parallel RepeatMasker runs.

    Sequences are tiled into non-overlapping regions of approximately equal total bp.
    Large scaffolds/chromosomes are split at chunk boundaries.

    Each batch produces:
      - A FASTA file with the extracted sequences
      - A BED file recording the original genomic coordinates (used downstream
        to restore global coordinates after RepeatMasker)

    COORDINATE-NAME CONTRACT:
      seqkit subseq --bed matches the BED chrom column against the sequence *ID*
      (--id-regexp default '^(\\S+)\\s?'), so the names written into the BED and
      into ${assembly_id}_seq_lengths.tsv must be IDs, not full headers.
      This is why fx2tab is called with -i/--only-id. Without it, an assembly
      with descriptions in its headers produces BED chrom names containing
      spaces, seqkit matches nothing, exits 0, and every batch FASTA is silently
      empty. build_coord_map.py and MERGE_RM_RESULTS also key their scaffold
      length lookups on column 1 of this file and require IDs.

    Inputs:
        - assembly_id:   Sample identifier
        - assembly_path: Path to genome assembly FASTA (ID-normalized upstream)
        - n_chunks:      Number of chunks to create (val, from params)

    Outputs:
        - batches:      path("batches/batch_*.fa")   — one FASTA per chunk
        - beds:         path("beds/batch_*.bed")     — matching BED per chunk
        - seq_lengths:  ${assembly_id}_seq_lengths.tsv (ID, length)

    Tools: seqkit (conda: bioconda::seqkit)
========================================================================================
*/

process PREP_RM_BATCHES {
    tag "${assembly_id}"
    label 'prep_rm_batches'

    input:
    tuple val(assembly_id), path(assembly_path)
    val n_chunks

    output:
    tuple val(assembly_id), path("batches/batch_*.fa"),  emit: batches
    tuple val(assembly_id), path("beds/batch_*.bed"),    emit: beds
    tuple val(assembly_id), path("${assembly_id}_seq_lengths.tsv"),     emit: seq_lengths

    script:
    """
    set -euo pipefail

    mkdir -p batches beds

    # ── 0. Assert the ID contract established by FILTER_ASSEMBLY ──────────────
    if grep -q '^>.*[[:space:]]' "${assembly_path}"; then
        echo "ERROR: ${assembly_path} contains FASTA headers with whitespace." >&2
        echo "       Sequence names must be bare IDs; run this assembly through" >&2
        echo "       FILTER_ASSEMBLY (header normalization) first." >&2
        grep -m1 '^>.*[[:space:]]' "${assembly_path}" >&2
        exit 1
    fi

    # ── 1. Get scaffold IDs and lengths ───────────────────────────────────────
    #   -i/--only-id is REQUIRED here; -n alone prints the full header.
    seqkit fx2tab -n -i -l "${assembly_path}" > ${assembly_id}_seq_lengths.tsv

    # ── 2. Generate per-chunk BED files ───────────────────────────────────────
    #   Walk scaffolds in input order, filling each chunk to ~target_size bp.
    #   Splits within scaffolds when a scaffold is larger than the remaining
    #   space in the current chunk. The last chunk absorbs any remainder.
    awk -v n="${n_chunks}" '
    BEGIN { FS = OFS = "\\t" }

    # Read scaffold names and lengths
    {
        names[NR] = \$1
        lens[NR]  = int(\$2)
        total    += lens[NR]
        n_seqs    = NR
    }

    END {
        chunk_size = int(total / n)
        if (chunk_size < 1) chunk_size = 1

        chunk  = 0
        filled = 0
        pos    = 0       # position within current scaffold
        si     = 1       # scaffold index

        while (si <= n_seqs) {
            remain_seq   = lens[si] - pos
            remain_chunk = chunk_size - filled

            # Last chunk: take everything that is left
            if (chunk == n - 1) remain_chunk = remain_seq

            if (remain_seq <= remain_chunk) {
                # Rest of scaffold fits in current chunk
                fname = sprintf("beds/batch_%04d.bed", chunk)
                print names[si], pos, lens[si] > fname
                filled += remain_seq
                si++
                pos = 0
            } else {
                # Split scaffold at chunk boundary
                end = pos + remain_chunk
                fname = sprintf("beds/batch_%04d.bed", chunk)
                print names[si], pos, end > fname
                pos    = end
                filled = 0
                chunk++
            }

            # Advance to next chunk when full (unless last)
            if (filled >= chunk_size && chunk < n - 1) {
                filled = 0
                chunk++
            }
        }
    }
    ' ${assembly_id}_seq_lengths.tsv

    # ── 3. Extract batch FASTAs, asserting a 1:1 BED-record -> sequence map ───
    for bed in beds/batch_*.bed; do
        batch_id=\$(basename "\$bed" .bed)
        seqkit subseq --bed "\$bed" "${assembly_path}" > "batches/\${batch_id}.fa"

        n_expected=\$(grep -c . "\$bed" || true)
        n_extracted=\$(grep -c '^>' "batches/\${batch_id}.fa" || true)

        if [ "\${n_extracted}" -ne "\${n_expected}" ]; then
            echo "ERROR: \${batch_id}: BED requested \${n_expected} regions but seqkit extracted \${n_extracted}." >&2
            echo "       BED chrom names must match FASTA sequence IDs exactly." >&2
            echo "       First BED record:" >&2
            head -1 "\$bed" >&2
            echo "       First FASTA header:" >&2
            grep -m1 '^>' "${assembly_path}" >&2
            exit 1
        fi
    done

    echo "PREP_RM_BATCHES ${assembly_id}: \$(ls -1 batches/batch_*.fa | wc -l) batches, all non-empty" >&2
    """

    stub:
    """
    mkdir -p batches beds
    touch ${assembly_id}_seq_lengths.tsv
    for i in \$(seq 0 \$((${n_chunks} - 1))); do
        printf -v pad "%04d" "\$i"
        touch "batches/batch_\${pad}.fa"
        touch "beds/batch_\${pad}.bed"
    done
    """
}
