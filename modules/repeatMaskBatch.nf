/*
========================================================================================
    REPEAT MASK BATCH
========================================================================================
    Run RepeatMasker on a single batch chunk, then adjust coordinates back to
    global genomic positions.

    seqkit subseq --bed outputs sequences in genome fasta order, NOT in BED
    file order. The coordinate mapping is therefore built by parsing the actual
    FASTA headers produced by seqkit (format: >scaffold_name_start-end:.)
    rather than assuming BED line N corresponds to FASTA sequence N.

    RepeatMasker enforces a 50-character limit on FASTA identifiers.
    Because the encoded names can exceed this, we:
      0. Rename sequences to short temp IDs and keep a mapping file
         (built from actual fasta header order, not BED line order)
      1. Run RepeatMasker on the short-named FASTA
      2. Map short IDs back to original scaffold + offset via header-derived map
      3. Adjust coordinates in .out, .masked, and .align files

    Inputs:
        - assembly_id:      Sample identifier
        - batch_id:         Batch identifier (e.g. "batch_0000")
        - batch_fasta:      FASTA file for this batch
        - batch_bed:        BED file with original genomic coordinates
        - combined_lib:     Combined repeat library for THIS assembly
        - seq_lengths:      Scaffold lengths TSV for THIS assembly

    Outputs:
        - adjusted .out, .fa, and .align files with global coordinates

    Tools: RepeatMasker (dfam/tetools container)
========================================================================================
*/

process REPEAT_MASK_BATCH {
    tag "${assembly_id}:${batch_id}"
    label 'repeat_mask_batch'
    //publishDir "${params.outdir}/repeats/chunks", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), val(batch_id), path(batch_fasta), path(batch_bed),
          path(combined_lib), path(seq_lengths)
    path(build_coord_map_script)
    path(verify_coord_map_script)

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
    awk '
    /^>/ { short = sprintf("seq_%04d", n++); print ">" short; next }
    { print }
    ' ${batch_fasta} > short_ids.fa

    # ── Build coordinate mapping from FASTA headers ───────────────────────────
    #   seqkit subseq --bed header format: >scaffold_name_start-end:.
    #   We parse these directly to get the correct seq_N → scaffold mapping,
    #   since seqkit outputs in genome order not BED order.
    python3 ${build_coord_map_script} ${batch_fasta} ${seq_lengths} coord_map.tsv

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
            map_offset[\$1] = int(\$3) - 1
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
            if (slen > 0) { \$8 = "(" slen - \$7 ")" }
        }
        print
    }
    ' short_ids.fa.out > ${assembly_id}.${batch_id}.adjusted.out

    # ── 3. Restore original names in masked FASTA ─────────────────────────────
    #   Headers become: >scaffold_name__start1__end
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

    # ── 4. Adjust .align file ─────────────────────────────────────────────────
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
    /^[[:space:]]*[0-9]/ && NF >= 14 {
        sid = \$5
        if (sid in map_scaff) {
            \$5 = map_scaff[sid]
            offset = map_offset[sid]
            slen   = map_slen[sid]
            \$6 = \$6 + offset
            \$7 = \$7 + offset
            if (slen > 0) { \$8 = "(" slen - \$7 ")" }
        }
        print
        next
    }
    { print }
    ' short_ids.fa.align > ${assembly_id}.${batch_id}.adjusted.align

    # ── 5. Verify coord_map correctness ──────────────────────────────────────
    python3 ${verify_coord_map_script} ${assembly_id}.${batch_id}.adjusted.fa >&2
    """

    stub:
    """
    touch ${assembly_id}.${batch_id}.adjusted.out
    touch ${assembly_id}.${batch_id}.adjusted.fa
    touch ${assembly_id}.${batch_id}.adjusted.align
    """
}