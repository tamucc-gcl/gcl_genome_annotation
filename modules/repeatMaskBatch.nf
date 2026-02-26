/*
========================================================================================
    REPEAT MASK BATCH
========================================================================================
    Run RepeatMasker on a single batch chunk, then adjust coordinates back to
    global genomic positions.

    seqkit subseq --bed names sequences as:
        scaffold_name:start-end
    which encodes the original scaffold name, 1-based start, and end separated
    by colons and dashes.

    RepeatMasker enforces a 50-character limit on FASTA identifiers.
    Because the encoded names can exceed this, we:
      0. Rename sequences to short temp IDs (seq_0000, seq_0001, ...)
         and keep a mapping file for restoration
      1. Run RepeatMasker on the short-named FASTA
      2. Restore original names in .out and .masked files
      3. Parse encoded names to extract original scaffold name and offset
      4. Add the offset to begin/end positions in the .out file
      5. Recalculate the (left) column using full scaffold lengths
      6. Restore original scaffold names in both .out and masked FASTA

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
    publishDir "${params.outdir}/repeats/chunks", mode: params.publish_dir_mode

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

    # ── 0. Rename sequences to short IDs (RepeatMasker max ID = 50 chars) ─────
    awk '
    /^>/ {
        orig = substr(\$0, 2)
        short = sprintf("seq_%04d", n++)
        map[short] = orig
        print ">" short
        next
    }
    { print }
    END {
        for (s in map) print s "\\t" map[s] > "id_map.tsv"
    }
    ' ${batch_fasta} > short_ids.fa

    # ── 1. Run RepeatMasker ───────────────────────────────────────────────────
    RepeatMasker \\
        -parallel ${rm_threads} \\
        -engine ncbi \\
        ${speed_flag} \\
        -cutoff 225 \\
        -lib ${combined_lib} \\
        -xsmall -small \\
        -gff \\
        short_ids.fa

    # ── 2. Restore original names in .out, then adjust to global coords ───────
    #
    # seqkit subseq --bed produces headers like:
    #   scaffold_1_sliding:7785837-15571672:.
    # Format: {scaffold_name}:{1-based_start}-{end}  (with optional trailing ":.")
    #
    # We load id_map.tsv and seq_lengths in BEGIN via getline so that awk
    # only reads the .out file as its main input. This avoids multi-file
    # NR/FNR portability issues across mawk/gawk.
    #
    awk -v mapfile="id_map.tsv" -v lenfile="${seq_lengths}" '
    BEGIN {
        while ((getline line < mapfile) > 0) {
            split(line, a, "\\t")
            id_map[a[1]] = a[2]
        }
        close(mapfile)

        while ((getline line < lenfile) > 0) {
            split(line, a, "\\t")
            scaff_len[a[1]] = int(a[2])
        }
        close(lenfile)
    }

    # Header lines (first 3 lines of RepeatMasker .out)
    NR <= 3 { print; next }

    {
        # Restore original long name from short ID
        if (\$5 in id_map) \$5 = id_map[\$5]

        seqname = \$5

        # Strip trailing ":." if present (seqkit subseq suffix)
        sub(/:\\.\$/, "", seqname)

        # Parse: scaffold_name:start-end
        # Split on ":" — last field contains "start-end"
        n_colons = split(seqname, colon_parts, ":")

        if (n_colons >= 2) {
            coord_field = colon_parts[n_colons]

            # Reconstruct original scaffold name from everything before last ":"
            orig_name = colon_parts[1]
            for (i = 2; i < n_colons; i++) {
                orig_name = orig_name ":" colon_parts[i]
            }

            # Parse start-end
            split(coord_field, se, "-")
            offset = int(se[1]) - 1
        } else {
            orig_name = seqname
            offset = 0
        }

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
    ' short_ids.fa.out > ${assembly_id}.${batch_id}.adjusted.out

    # ── 3. Restore original names in masked FASTA, then adjust headers ────────
    #
    # Convert: >scaffold_name:start-end:.  →  >scaffold_name__start__end
    #
    awk '
    NR == FNR { id_map[\$1] = \$2; next }
    /^>/ {
        short = substr(\$1, 2)
        header = (short in id_map) ? id_map[short] : short

        # Strip trailing ":."
        sub(/:\\.\$/, "", header)

        # Parse scaffold_name:start-end
        n_colons = split(header, colon_parts, ":")

        if (n_colons >= 2) {
            coord_field = colon_parts[n_colons]
            orig_name = colon_parts[1]
            for (i = 2; i < n_colons; i++) {
                orig_name = orig_name ":" colon_parts[i]
            }
            split(coord_field, se, "-")
            print ">" orig_name "__" int(se[1]) "__" int(se[2])
        } else {
            print ">" header "__1__1"
        }
        next
    }
    { print }
    ' id_map.tsv short_ids.fa.masked > ${assembly_id}.${batch_id}.adjusted.fa
    """

    stub:
    """
    touch ${assembly_id}.${batch_id}.adjusted.out
    touch ${assembly_id}.${batch_id}.adjusted.fa
    """
}