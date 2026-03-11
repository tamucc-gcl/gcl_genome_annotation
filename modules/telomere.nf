/*
========================================================================================
    TELOMERE ANNOTATION MODULE
========================================================================================
    Identify telomeric repeat regions using tidk (Telomere Identification toolKit)
    and convert windowed counts to GFF3 annotation.

    tidk search scans for a user-specified repeat motif in windows across the
    genome and reports per-window counts. This module:
      1. Runs tidk search for forward/reverse motif counts
      2. Generates an SVG plot via tidk plot for QC
      3. Converts windowed counts to discrete GFF3 regions by merging
         consecutive above-threshold windows and classifying by position
         (5_prime, 3_prime, interstitial)

    Inputs:
        - assembly_id:      Sample identifier
        - assembly_fasta:   Genome assembly FASTA (unmasked preferred)

    Outputs:
        - telomere_gff3:    ${assembly_id}.telomeres.gff3   (GFF3 annotation)
        - telomere_tsv:     ${assembly_id}.telomeres.tsv    (tidk raw windowed counts)
        - telomere_plot:    ${assembly_id}.telomeres.svg    (tidk density plot)

    Tools: tidk (conda: bioconda::tidk), Python 3 (stdlib only)

    Reference:
        Brown MR, Gonzalez de La Rosa PM, Blaxter M. tidk: a toolkit to rapidly
        identify telomeric repeats from genomic datasets. Bioinformatics. 2025;
        41(2):btaf049. doi:10.1093/bioinformatics/btaf049
========================================================================================
*/

process TELOMERE_ANNOTATION {
    tag "${assembly_id}"
    label 'telomere'
    publishDir "${params.outdir}/telomeres", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(assembly_fasta)

    output:
    tuple val(assembly_id), path("${assembly_id}.telomeres.gff3"),  emit: telomere_gff3
    tuple val(assembly_id), path("${assembly_id}.telomeres.tsv"),   emit: telomere_tsv
    tuple val(assembly_id), path("${assembly_id}.telomeres.svg"),   emit: telomere_plot

    script:
    def motif       = params.telomere_motif ?: 'TTAGGG'
    def window_size = params.telomere_window ?: 10000
    def min_count   = params.telomere_min_count ?: 2
    def end_window  = params.telomere_end_window ?: 50000
    """
    set -euo pipefail

    # ── 1. Run tidk search — windowed telomere repeat counts ──────────────────
    mkdir -p tidk_out

    tidk search \\
        --string ${motif} \\
        --window ${window_size} \\
        --output ${assembly_id} \\
        --dir tidk_out \\
        ${assembly_fasta}

    cp tidk_out/${assembly_id}_telomeric_repeat_windows.tsv ${assembly_id}.telomeres.tsv

    # ── 2. Generate SVG plot for QC ───────────────────────────────────────────
    tidk plot \\
        --tsv ${assembly_id}.telomeres.tsv \\
        --output ${assembly_id}.telomeres

    # ── 3. Convert windowed counts to GFF3 ────────────────────────────────────
    python3 - <<'PYEOF'
import sys

# ── Parameters ────────────────────────────────────────────────────────────
WINDOW_SIZE = ${window_size}
MIN_COUNT   = ${min_count}   # min total repeats (fwd + rev) in a window to call it telomeric
END_WINDOW  = ${end_window}  # bp from scaffold end to classify as terminal
MOTIF       = "${motif}"
ASM         = "${assembly_id}"

# ── Read tidk TSV ─────────────────────────────────────────────────────────
#   Columns: id  window  forward_repeat_number  reverse_repeat_number  telomeric_repeat
scaffold_windows = {}   # scaffold → [(window_start, fwd, rev), ...]
scaffold_lengths = {}   # scaffold → max window_start + WINDOW_SIZE

with open(f"{ASM}.telomeres.tsv") as fh:
    header = fh.readline()
    for line in fh:
        line = line.strip()
        if not line:
            continue
        cols = line.split("\\t")
        scaffold = cols[0]
        win_start = int(cols[1])
        fwd = int(cols[2])
        rev = int(cols[3])

        scaffold_windows.setdefault(scaffold, []).append((win_start, fwd, rev))

        # Track max coordinate to estimate scaffold length
        end = win_start + WINDOW_SIZE
        if scaffold not in scaffold_lengths or end > scaffold_lengths[scaffold]:
            scaffold_lengths[scaffold] = end

# ── Identify telomeric regions ────────────────────────────────────────────
def classify_position(start, end, seq_len, window):
    at_start = start <= window
    at_end   = end >= seq_len - window
    if at_start and at_end:
        return "spanning"
    if at_start:
        return "5_prime"
    if at_end:
        return "3_prime"
    return "interstitial"

gff_lines = []
telo_counter = 0

for scaffold in sorted(scaffold_windows.keys()):
    windows = sorted(scaffold_windows[scaffold], key=lambda x: x[0])
    seq_len = scaffold_lengths[scaffold]

    # Find runs of consecutive windows above threshold
    in_region = False
    region_start = 0
    region_fwd = 0
    region_rev = 0
    prev_end = 0

    for win_start, fwd, rev in windows:
        total = fwd + rev
        win_end = win_start + WINDOW_SIZE

        if total >= MIN_COUNT:
            if not in_region:
                # Start new region
                in_region = True
                region_start = win_start
                region_fwd = fwd
                region_rev = rev
            else:
                # Extend region
                region_fwd += fwd
                region_rev += rev
            prev_end = win_end
        else:
            if in_region:
                # Close region
                telo_counter += 1
                region_end = prev_end
                n_repeats = region_fwd + region_rev
                strand = "+" if region_fwd >= region_rev else "-"
                position = classify_position(
                    region_start, region_end, seq_len, END_WINDOW
                )

                feat_id = f"{ASM}_telo{telo_counter}"
                attrs = (
                    f"ID={feat_id};Name=telomeric_repeat;"
                    f"motif={MOTIF};n_repeats={n_repeats};"
                    f"forward_count={region_fwd};reverse_count={region_rev};"
                    f"position={position}"
                )

                # GFF3: score = total repeat count
                gff_lines.append(
                    f"{scaffold}\\ttidk\\ttelomeric_repeat\\t"
                    f"{region_start + 1}\\t{region_end}\\t"
                    f"{n_repeats}\\t{strand}\\t.\\t{attrs}"
                )

                in_region = False
                region_fwd = 0
                region_rev = 0

    # Close final region if scaffold ends in telomere
    if in_region:
        telo_counter += 1
        region_end = prev_end
        n_repeats = region_fwd + region_rev
        strand = "+" if region_fwd >= region_rev else "-"
        position = classify_position(
            region_start, region_end, seq_len, END_WINDOW
        )

        feat_id = f"{ASM}_telo{telo_counter}"
        attrs = (
            f"ID={feat_id};Name=telomeric_repeat;"
            f"motif={MOTIF};n_repeats={n_repeats};"
            f"forward_count={region_fwd};reverse_count={region_rev};"
            f"position={position}"
        )

        gff_lines.append(
            f"{scaffold}\\ttidk\\ttelomeric_repeat\\t"
            f"{region_start + 1}\\t{region_end}\\t"
            f"{n_repeats}\\t{strand}\\t.\\t{attrs}"
        )

# ── Write GFF3 ────────────────────────────────────────────────────────────
with open(f"{ASM}.telomeres.gff3", "w") as out:
    out.write("##gff-version 3\\n")
    out.write(f"# Telomere annotation for {ASM} (tidk search)\\n")
    out.write(f"# Motif: {MOTIF}, window: {WINDOW_SIZE} bp, min_count: {MIN_COUNT}\\n")
    out.write(f"# Total telomeric regions: {telo_counter}\\n")
    for line in gff_lines:
        out.write(line + "\\n")

# ── Stats ─────────────────────────────────────────────────────────────────
n_5p  = sum(1 for l in gff_lines if "position=5_prime" in l)
n_3p  = sum(1 for l in gff_lines if "position=3_prime" in l)
n_int = sum(1 for l in gff_lines if "position=interstitial" in l)
n_span = sum(1 for l in gff_lines if "position=spanning" in l)
print(f"Telomere annotation: {telo_counter} regions "
      f"(5'={n_5p}, 3'={n_3p}, interstitial={n_int}, spanning={n_span})",
      file=sys.stderr)
PYEOF
    """

    stub:
    """
    touch ${assembly_id}.telomeres.gff3
    touch ${assembly_id}.telomeres.tsv
    touch ${assembly_id}.telomeres.svg
    """
}