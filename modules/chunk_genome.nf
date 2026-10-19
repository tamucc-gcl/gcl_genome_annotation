/*
========================================================================================
    Chunk Genome MODULE
========================================================================================
    Split masked genome into N balanced chunks for parallel processing (e.g. Infernal).

    Uses greedy bin-packing (longest-first) to assign whole scaffolds to bins
    so that total bp per chunk is approximately equal. This avoids the load
    imbalance of naive sequential splitting when scaffold sizes vary widely
    (e.g. chromosome-scale assemblies).

    Algorithm:
      1. Get scaffold names + lengths
      2. Sort scaffolds by length descending
      3. Assign each scaffold to the bin with the smallest current total
      4. Extract per-bin FASTAs using seqkit grep

    Inputs:
        - assembly_id:      Sample identifier
        - masked_assembly:  Soft-masked genome assembly FASTA

    Outputs:
        - chunks:  tuple val(assembly_id), path("${assembly_id}_genome_chunks/*")

    Tools: seqkit (conda: bioconda::seqkit)
========================================================================================
*/

process CHUNK_GENOME {
    tag "${assembly_id}"
    label 'chunk_genome'

    input:
    tuple val(assembly_id), path(masked_assembly)

    output:
    tuple val(assembly_id), path("${assembly_id}_genome_chunks/chunk_*.fa"), emit: chunks
    tuple val(assembly_id), path("seq_lengths.tsv"), emit: lengths

    script:
    """
    set -euo pipefail

    mkdir -p ${assembly_id}_genome_chunks

    # ── 1. Get scaffold names and lengths ─────────────────────────────────────
    seqkit fx2tab -n -l "${masked_assembly}" > seq_lengths.tsv

    # ── 2. Greedy bin-packing: assign scaffolds to N balanced chunks ──────────
    #   Sort by length desc, assign each to the bin with smallest current total.
    #   Outputs one file per bin: chunk_000.names, chunk_001.names, ...
    awk -v n="${params.n_infernal_chunks}" '
    BEGIN {
        for (i = 0; i < n; i++) bin_total[i] = 0
    }
    {
        names[NR] = \$1
        lens[NR]  = int(\$2)
        count = NR
    }
    END {
        # Sort indices by length descending (simple selection sort — fine for <1000 scaffolds)
        for (i = 1; i <= count; i++) order[i] = i
        for (i = 1; i <= count; i++) {
            max_idx = i
            for (j = i + 1; j <= count; j++) {
                if (lens[order[j]] > lens[order[max_idx]]) max_idx = j
            }
            if (max_idx != i) {
                tmp = order[i]; order[i] = order[max_idx]; order[max_idx] = tmp
            }
        }

        # Greedy assignment: put each scaffold in the lightest bin
        for (k = 1; k <= count; k++) {
            idx = order[k]

            # Find bin with smallest total
            min_bin = 0
            for (b = 1; b < n; b++) {
                if (bin_total[b] < bin_total[min_bin]) min_bin = b
            }

            bin_total[min_bin] += lens[idx]
            fname = sprintf("chunk_%03d.names", min_bin)
            print names[idx] > fname
        }

        # Print bin sizes for logging
        for (i = 0; i < n; i++) {
            printf "Chunk %03d: %d bp\\n", i, bin_total[i] > "/dev/stderr"
        }
    }
    ' seq_lengths.tsv

    # ── 3. Extract per-chunk FASTAs ───────────────────────────────────────────
    for namefile in chunk_*.names; do
        chunk_id=\${namefile%.names}
        seqkit grep -f "\$namefile" "${masked_assembly}" \
            > "${assembly_id}_genome_chunks/\${chunk_id}.fa"
    done
    """

    stub:
    """
    mkdir -p ${assembly_id}_genome_chunks
    touch seq_lengths.tsv
    for i in \$(seq 0 \$((${params.n_infernal_chunks} - 1))); do
        printf -v pad "%03d" "\$i"
        touch "${assembly_id}_genome_chunks/chunk_\${pad}.fa"
    done
    """
}