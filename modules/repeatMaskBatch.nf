/*
========================================================================================
    REPEAT MASK BATCH
========================================================================================
    Run RepeatMasker on a single batch chunk, then adjust coordinates back to
    global genomic positions.

    seqkit subseq --bed names sequences as:
        scaffold_name_start-end:.
    which encodes: {original_name}_{1-based_start}-{end}:.

    The trailing ":." is a seqkit artifact. The coordinate token is the last
    underscore-delimited field that contains a dash (e.g. "7785837-15571672").

    RepeatMasker enforces a 50-character limit on FASTA identifiers.
    Because the encoded names can exceed this, we:
      0. Rename sequences to short temp IDs and keep a mapping file
      1. Run RepeatMasker on the short-named FASTA
      2. Restore original names, parse encoded coords, adjust to global positions

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

    # ── 2. Build a direct lookup: short_id → scaffold_name, offset ────────────
    #
    # id_map.tsv maps short IDs to seqkit-encoded names like:
    #   seq_0000  scaffold_1_sliding_1-7785836:.
    #
    # We parse these to extract the original scaffold name and the 1-based
    # start coordinate (offset). The coordinate token is the last
    # underscore-delimited field that contains a "-".
    #
    # Output: parsed_map.tsv with columns: short_id, scaffold_name, offset
    #
    awk -v lenfile="${seq_lengths}" '
    BEGIN {
        # Load scaffold lengths
        while ((getline line < lenfile) > 0) {
            split(line, a, "\\t")
            scaff_len[a[1]] = int(a[2])
        }
        close(lenfile)
    }
    {
        short = \$1
        encoded = \$2

        # Strip trailing ":." if present
        sub(/:\\.\$/, "", encoded)

        # Split on "_" to find the coordinate token
        n = split(encoded, parts, "_")

        # Find the last part that contains a "-" — that is the coord token
        coord_idx = 0
        for (i = n; i >= 1; i--) {
            if (index(parts[i], "-") > 0) {
                coord_idx = i
                break
            }
        }

        if (coord_idx > 0) {
            # Reconstruct scaffold name from parts before the coord token
            orig_name = parts[1]
            for (i = 2; i < coord_idx; i++) {
                orig_name = orig_name "_" parts[i]
            }

            # Parse start-end
            split(parts[coord_idx], se, "-")
            offset = int(se[1]) - 1
        } else {
            # No coordinate encoding found
            orig_name = encoded
            offset = 0
        }

        # Store for lookup
        print short, orig_name, offset, (orig_name in scaff_len ? scaff_len[orig_name] : 0)
    }
    ' id_map.tsv > parsed_map.tsv

    # ── 3. Adjust .out file using parsed mapping ──────────────────────────────
    awk '
    BEGIN {
        while ((getline line < "parsed_map.tsv") > 0) {
            split(line, a, " ")
            short_id = a[1]
            map_name[short_id]   = a[2]
            map_offset[short_id] = int(a[3])
            map_scafflen[short_id] = int(a[4])
        }
        close("parsed_map.tsv")
    }

    # Header lines
    NR <= 3 { print; next }

    {
        sid = \$5

        if (sid in map_name) {
            \$5 = map_name[sid]
            offset = map_offset[sid]
            slen   = map_scafflen[sid]

            \$6 = \$6 + offset
            \$7 = \$7 + offset

            if (slen > 0) {
                \$8 = "(" slen - \$7 ")"
            }
        }

        print
    }
    ' short_ids.fa.out > ${assembly_id}.${batch_id}.adjusted.out

    # ── 4. Restore original names in masked FASTA ─────────────────────────────
    #
    # Convert encoded headers to: >scaffold_name__start__end
    # (double-underscore delimited for MERGE_RM_RESULTS to parse)
    #
    awk '
    NR == FNR {
        encoded = \$2
        sub(/:\\.\$/, "", encoded)

        n = split(encoded, parts, "_")
        coord_idx = 0
        for (i = n; i >= 1; i--) {
            if (index(parts[i], "-") > 0) { coord_idx = i; break }
        }

        if (coord_idx > 0) {
            orig = parts[1]
            for (i = 2; i < coord_idx; i++) orig = orig "_" parts[i]
            split(parts[coord_idx], se, "-")
            fa_name[\$1] = orig "__" int(se[1]) "__" int(se[2])
        } else {
            fa_name[\$1] = encoded "__1__1"
        }
        next
    }
    /^>/ {
        short = substr(\$1, 2)
        if (short in fa_name) print ">" fa_name[short]
        else print
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