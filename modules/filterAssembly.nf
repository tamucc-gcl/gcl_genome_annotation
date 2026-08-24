/*
========================================================================================
    FILTER ASSEMBLY
========================================================================================
    Single point of entry for every assembly in the pipeline.

    Responsibilities (in order):
      1. NORMALIZE HEADERS. Every downstream process — PREP_RM_BATCHES,
         build_coord_map.py, MERGE_RM_RESULTS, BRAKER, AGAT, every GFF3 seqid —
         assumes the sequence name is the first whitespace-delimited token and
         contains no characters that are special in GFF3/BED/SAM. NCBI
         `_genomic.fna` downloads violate this (`>ACC.1 Genus species isolate ...`),
         as do many assembler outputs (`>ctg1|arrow|pilon`). Normalizing here,
         once, means no downstream module ever has to think about it.
      2. Filter contigs by minimum length.
      3. Emit an auditable clean_id -> original_header map so the normalization
         is reversible and traceable.

    CONTRACT ESTABLISHED BY THIS PROCESS:
        Sequence names in ${assembly_id}.filtered.fasta match ^[A-Za-z0-9._-]+$
        and are unique. Everything downstream may rely on this.

    Inputs:
        - assembly_id:       Sample identifier
        - assembly_path:     Path to genome assembly FASTA (any header style)
        - min_contig_length: Minimum contig length to retain (params)

    Outputs:
        - filtered_assembly: ${assembly_id}.filtered.fasta
        - header_map:        ${assembly_id}.header_map.tsv
                             (clean_id, original_id, original_header, status)
        - assembly_stats:    ${assembly_id}.assembly_stats.tsv

    Tools: seqkit (conda: bioconda::seqkit), assembly-stats (conda: bioconda::assembly-stats)
========================================================================================
*/

process FILTER_ASSEMBLY {
    tag "${assembly_id}"
    label 'filter_assembly'
    publishDir "${params.outdir}/filtered_assemblies", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(assembly_path)

    output:
    tuple val(assembly_id), path("${assembly_id}.filtered.fasta"),     emit: filtered_assembly
    tuple val(assembly_id), path("${assembly_id}.header_map.tsv"),     emit: header_map
    tuple val(assembly_id), path("${assembly_id}.assembly_stats.tsv"), emit: assembly_stats

    script:
    """
    set -euo pipefail

    # ── 1. Normalize headers ──────────────────────────────────────────────────
    #   >JBKOQA010000001.1 Taeniamia zosterophora isolate ... shotgun sequence
    #   becomes
    #   >JBKOQA010000001.1
    #
    #   Anything outside [A-Za-z0-9._-] in the ID itself becomes '_'.
    #   Collisions are a hard error: silently merging two scaffolds into one
    #   name would corrupt every coordinate downstream.
    awk '
    /^>/ {
        hdr = substr(\$0, 2)
        sub(/^[[:space:]]+/, "", hdr)

        orig_id = hdr
        sub(/[[:space:]].*\$/, "", orig_id)

        if (orig_id == "") {
            print "ERROR: empty sequence ID at input line " NR > "/dev/stderr"
            exit 1
        }

        clean = orig_id
        gsub(/[^A-Za-z0-9._-]/, "_", clean)

        if (clean in seen) {
            print "ERROR: sequence ID collision after normalization: [" clean "]" > "/dev/stderr"
            print "       produced by both [" seen[clean] "] and [" orig_id "]" > "/dev/stderr"
            exit 1
        }
        seen[clean] = orig_id

        if (length(clean) > 50) {
            print "WARNING: sequence ID exceeds 50 characters (RepeatMasker limit): [" clean "]" > "/dev/stderr"
        }

        print clean "\\t" orig_id "\\t" hdr > "header_map.raw.tsv"
        print ">" clean
        next
    }
    { print }
    ' "${assembly_path}" > renamed.fasta

    # ── 2. Filter contigs by minimum length ───────────────────────────────────
    #   NOTE: -g/--remove-gaps is deliberately NOT used. Its default gap letters
    #   are " -." so it silently deletes those characters from the sequence,
    #   shifting every downstream coordinate.
    seqkit seq -m ${params.min_contig_length} renamed.fasta > ${assembly_id}.filtered.fasta

    # ── 3. Fail loudly rather than emitting an empty FASTA ────────────────────
    total_before=\$(grep -c '^>' renamed.fasta || true)
    total_after=\$(grep -c '^>' ${assembly_id}.filtered.fasta || true)

    if [ "\${total_before}" -eq 0 ]; then
        echo "ERROR: no sequences found in ${assembly_path}" >&2
        exit 1
    fi
    if [ "\${total_after}" -eq 0 ]; then
        echo "ERROR: no contigs >= ${params.min_contig_length} bp survived filtering (input had \${total_before})" >&2
        exit 1
    fi

    removed=\$((total_before - total_after))
    echo "Filtered ${assembly_id}: kept \${total_after}/\${total_before} contigs (removed \${removed} < ${params.min_contig_length} bp)" >&2

    # ── 4. Header map with retention status ───────────────────────────────────
    seqkit fx2tab -n -i ${assembly_id}.filtered.fasta | cut -f1 | sort -u > kept_ids.txt

    awk -F'\\t' '
    BEGIN { OFS = "\\t"; print "clean_id", "original_id", "original_header", "status" }
    NR == FNR { kept[\$1] = 1; next }
    { print \$1, \$2, \$3, (\$1 in kept ? "retained" : "removed") }
    ' kept_ids.txt header_map.raw.tsv > ${assembly_id}.header_map.tsv

    # ── 5. Generate assembly stats ────────────────────────────────────────────
    assembly-stats -t ${assembly_id}.filtered.fasta > ${assembly_id}.assembly_stats.tsv
    """

    stub:
    """
    touch ${assembly_id}.filtered.fasta
    touch ${assembly_id}.header_map.tsv
    touch ${assembly_id}.assembly_stats.tsv
    """
}
