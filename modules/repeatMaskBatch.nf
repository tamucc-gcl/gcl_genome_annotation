/*
========================================================================================
    REPEAT MASK BATCH
========================================================================================
    Run RepeatMasker on a single batch chunk, then adjust coordinates back to
    global genomic positions.

    The BED file from PREP_RM_BATCHES is the authoritative source for the
    mapping between batch sequences and original genomic coordinates.
    Scaffold names in the input assembly may contain colons, dashes, and
    underscores, making header parsing unreliable. Instead we use the BED
    file order (which matches the FASTA order from seqkit subseq --bed)
    to build a positional mapping: sequence N in the batch → BED line N.

    RepeatMasker enforces a 50-character limit on FASTA identifiers.
    Because the encoded names can exceed this, we:
      0. Rename sequences to short temp IDs and keep a mapping file
      1. Run RepeatMasker on the short-named FASTA
      2. Map short IDs back to original scaffold + offset via BED
      3. Adjust coordinates in .out and .masked files

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
        path("${assembly_id}.${batch_id}.adjusted.align"),
        emit: batch_results

    script:
    def speed_flag = params.job_speed ? "-${params.job_speed}" : "-s"
    def rm_threads = Math.max(1, (task.cpus / 4) as int)
    """
    set -euo pipefail

    # ── 0. Rename sequences to short IDs (RepeatMasker max ID = 50 chars) ─────
    #   Also record the order (seq_0000, seq_0001, ...) which matches BED line order
    awk '
    /^>/ {
        short = sprintf("seq_%04d", n++)
        print ">" short
        next
    }
    { print }
    ' ${batch_fasta} > short_ids.fa

    # ── Build coordinate mapping from BED ─────────────────────────────────────
    #   BED lines are in the same order as FASTA sequences produced by seqkit.
    #   BED format: scaffold_name  0-based_start  end
    #   We create: short_id  scaffold_name  1-based_start  end
    #
    #   Also load scaffold lengths for recalculating the (left) column.
    awk -v lenfile="${seq_lengths}" '
    BEGIN {
        while ((getline line < lenfile) > 0) {
            split(line, a, "\\t")
            scaff_len[a[1]] = int(a[2])
        }
        close(lenfile)
    }
    {
        short = sprintf("seq_%04d", NR - 1)
        scaff = \$1
        start1 = int(\$2) + 1   # Convert 0-based BED start to 1-based
        end_coord = int(\$3)
        slen = (scaff in scaff_len) ? scaff_len[scaff] : 0
        print short "\\t" scaff "\\t" start1 "\\t" end_coord "\\t" slen
    }
    ' ${batch_bed} > coord_map.tsv

    # ── 1. Run RepeatMasker ───────────────────────────────────────────────────
    RepeatMasker \\
        -parallel ${rm_threads} \\
        -engine ncbi \\
        ${speed_flag} \\
        -cutoff 225 \\
        -lib ${combined_lib} \\
        -xsmall -small \\
        -a \\
        -gff \\
        short_ids.fa

    # ── 2. Adjust .out file ───────────────────────────────────────────────────
    awk '
    BEGIN {
        FS = "\\t"
        while ((getline < "coord_map.tsv") > 0) {
            map_scaff[\$1]  = \$2
            map_offset[\$1] = int(\$3) - 1   # offset = start1 - 1
            map_end[\$1]    = int(\$4)
            map_slen[\$1]   = int(\$5)
        }
        close("coord_map.tsv")
        FS = " "
    }

    NR <= 3 { print; next }

    {
        sid = \$5

        if (sid in map_scaff) {
            \$5 = map_scaff[sid]
            offset = map_offset[sid]
            slen   = map_slen[sid]

            \$6 = \$6 + offset
            \$7 = \$7 + offset

            if (slen > 0) {
                \$8 = "(" slen - \$7 ")"
            }
        }

        print
    }
    ' short_ids.fa.out > ${assembly_id}.${batch_id}.adjusted.out

    # ── 3. Restore original names in masked FASTA ─────────────────────────────
    #   Headers become: >scaffold_name__start1__end
    #   (double-underscore delimited for MERGE_RM_RESULTS to parse)
    awk '
    BEGIN {
        while ((getline < "coord_map.tsv") > 0) {
            split(\$0, a, "\\t")
            fa_hdr[a[1]] = a[2] "__" a[3] "__" a[4]
        }
        close("coord_map.tsv")
    }
    /^>/ {
        short = substr(\$1, 2)
        if (short in fa_hdr) print ">" fa_hdr[short]
        else print
        next
    }
    { print }
    ' short_ids.fa.masked > ${assembly_id}.${batch_id}.adjusted.fa
    """

    # ── 5. Adjust .align file ─────────────────────────────────────────────────
    #   .align has summary lines that mirror .out format with scaffold names
    #   and coordinates. We replace short IDs and adjust coordinates on lines
    #   that start with a score (numeric first field) similar to .out format.
    #   Other lines (alignments, blank lines, headers) pass through unchanged.
    awk '
    BEGIN {
        while ((getline < "coord_map.tsv") > 0) {
            split(\$0, a, "\\t")
            map_scaff[a[1]]  = a[2]
            map_offset[a[1]] = int(a[3]) - 1
            map_slen[a[1]]   = int(a[5])
        }
        close("coord_map.tsv")
    }

    # Summary/annotation lines start with a number (SW score)
    /^[[:space:]]*[0-9]/ && NF >= 14 {
        sid = \$5
        if (sid in map_scaff) {
            \$5 = map_scaff[sid]
            offset = map_offset[sid]
            slen   = map_slen[sid]

            \$6 = \$6 + offset
            \$7 = \$7 + offset

            if (slen > 0) {
                \$8 = "(" slen - \$7 ")"
            }
        }
        print
        next
    }

    # All other lines pass through unchanged
    { print }
    ' short_ids.fa.align > ${assembly_id}.${batch_id}.adjusted.align
    """

    stub:
    """
    touch ${assembly_id}.${batch_id}.adjusted.out
    touch ${assembly_id}.${batch_id}.adjusted.fa
    touch ${assembly_id}.${batch_id}.adjusted.align
    """
}