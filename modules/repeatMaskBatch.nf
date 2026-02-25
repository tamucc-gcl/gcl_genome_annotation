/*
========================================================================================
    REPEAT MASK BATCH
========================================================================================
    Run RepeatMasker on a single batch chunk, then adjust coordinates back to
    global genomic positions.

    seqkit subseq --bed names sequences as:
        scaffold_1_84087028-93430037:.
    which encodes: {original_name}_{1-based_start}-{end}:.

    After RepeatMasker completes, this module:
      1. Parses encoded names to extract original scaffold name and offset
      2. Adds the offset to begin/end positions in the .out file
      3. Recalculates the (left) column using full scaffold lengths
      4. Restores original scaffold names in both .out and masked FASTA

    Inputs:
        - assembly_id:      Sample identifier
        - batch_id:         Batch identifier (e.g. "batch_0000")
        - batch_fasta:      FASTA file for this batch
        - batch_bed:        BED file with original genomic coordinates
        - combined_lib:     Combined repeat library for THIS assembly
        - seq_lengths:      Scaffold lengths TSV for THIS assembly

    Outputs:
        - adjusted .out and .fa files with global coordinates

    Tools: RepeatMasker (dfam/tetools container)
========================================================================================
*/

process REPEAT_MASK_BATCH {
    tag "${assembly_id}:${batch_id}"
    label 'repeat_mask_batch'

    input:
    tuple val(assembly_id), val(batch_id), path(batch_fasta), path(batch_bed),
          path(combined_lib), path(seq_lengths)

    output:
    tuple val(assembly_id), val(batch_id),
        path("${assembly_id}.${batch_id}.adjusted.out"),
        path("${assembly_id}.${batch_id}.adjusted.fa"),
        emit: batch_results

    script:
    def speed_flag = params.job_speed ? "-${params.job_speed}" : "-s"
    def rm_threads = Math.max(1, (task.cpus / 4) as int)
    """
    set -euo pipefail

    # ── 1. Run RepeatMasker ───────────────────────────────────────────────────
    RepeatMasker \\
        -parallel ${rm_threads} \\
        -engine ncbi \\
        ${speed_flag} \\
        -cutoff 225 \\
        -lib ${combined_lib} \\
        -xsmall -small \\
        -gff \\
        ${batch_fasta}

    # ── 2. Adjust .out coordinates to global ──────────────────────────────────
    awk '
    # First pass: load full scaffold lengths
    NR == FNR {
        scaff_len[\$1] = int(\$2)
        next
    }

    # Second pass: process the .out file
    FNR <= 3 { print; next }

    {
        seqname = \$5

        # Parse: last _DIGITS-DIGITS:. is the coordinate suffix
        n = split(seqname, parts, "_")
        tail = parts[n]

        # Reconstruct original name from all parts except the last
        orig_name = parts[1]
        for (i = 2; i < n; i++) {
            orig_name = orig_name "_" parts[i]
        }

        # Parse start from tail: "84087028-93430037:."
        split(tail, coords, "[-:]")
        offset = int(coords[1]) - 1

        # Adjust begin/end coordinates
        \$5 = orig_name
        \$6 = \$6 + offset
        \$7 = \$7 + offset

        # Recalculate (left) using full scaffold length
        if (orig_name in scaff_len) {
            \$8 = "(" scaff_len[orig_name] - \$7 ")"
        }

        print
    }
    ' ${seq_lengths} ${batch_fasta}.out > ${assembly_id}.${batch_id}.adjusted.out

    # ── 3. Adjust masked FASTA headers ────────────────────────────────────────
    awk '
    /^>/ {
        header = substr(\$0, 2)
        gsub(/:\\.\$/, "", header)

        n = split(header, parts, "_")
        tail = parts[n]

        orig_name = parts[1]
        for (i = 2; i < n; i++) {
            orig_name = orig_name "_" parts[i]
        }

        split(tail, coords, "-")
        print ">" orig_name "__" int(coords[1]) "__" int(coords[2])
        next
    }
    { print }
    ' ${batch_fasta}.masked > ${assembly_id}.${batch_id}.adjusted.fa
    """

    stub:
    """
    touch ${assembly_id}.${batch_id}.adjusted.out
    touch ${assembly_id}.${batch_id}.adjusted.fa
    """
}