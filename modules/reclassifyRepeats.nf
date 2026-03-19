/*
========================================================================================
    RECLASSIFY UNKNOWN REPEATS
========================================================================================
    Use Terrier (deep learning classifier trained on RepBase) to reclassify
    repeat families that RepeatClassifier labeled as "Unknown".

    RepeatClassifier relies on homology to known TE families in Dfam/RepBase.
    For non-model species with lineage-specific TEs, many families lack close
    database matches and are labeled "Unknown" even though they have clear TE
    structural signatures. Terrier uses k-mer composition and learned sequence
    patterns to classify these into Wicker superfamilies without requiring
    homology, reducing the Unknown fraction substantially.

    Already-classified families (those with a class/superfamily from
    RepeatClassifier) are trusted and passed through unchanged — homology-based
    classification is more precise when a close match exists.

    Strategy:
      1. Split classified library into Known (has class) and Unknown subsets
      2. Run Terrier on Unknown subset with user-defined confidence threshold
      3. Merge reclassified families back with Known families
      4. Produce a reclassification report

    Inputs:
        - assembly_id:       Sample identifier
        - classified_lib:    RepeatClassifier output FASTA (headers have #Class)

    Outputs:
        - reclassified_lib:  ${assembly_id}_reclassified_repeats.fa
        - reclass_report:    ${assembly_id}_reclassification_report.tsv

    Tools: Terrier (pip: terrier-bio), Python 3 (stdlib only), seqkit
    Reference:
        Turnbull et al. (2025) Terrier: a deep learning repeat classifier.
        Briefings in Bioinformatics 26(4):bbaf442
========================================================================================
*/

process RECLASSIFY_REPEATS {
    tag "${assembly_id}"
    label 'reclassify_repeats'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(classified_lib)

    output:
    tuple val(assembly_id), path("${assembly_id}_reclassified_repeats.fa"),
        emit: reclassified_lib
    tuple val(assembly_id), path("${assembly_id}_reclassification_report.tsv"),
        emit: reclass_report

    script:
    def threshold = params.terrier_threshold ?: 0.7
    """
    set -euo pipefail

    # ── 1. Split library into Known and Unknown ───────────────────────────────
    #   RepeatClassifier format: >famname#Class/Superfamily
    #   Unknown families have #Unknown or no classification after #
    python3 - <<'PYEOF'
import sys

known_f = open("known_repeats.fa", "w")
unknown_f = open("unknown_repeats.fa", "w")
unknown_ids = set()
known_ids = set()

current_header = None
current_seq = []
is_unknown = False

def flush():
    if current_header is None:
        return
    if is_unknown:
        unknown_f.write(current_header + "\\n")
        unknown_ids.add(current_header.split()[0][1:])
        for s in current_seq:
            unknown_f.write(s + "\\n")
    else:
        known_f.write(current_header + "\\n")
        known_ids.add(current_header.split()[0][1:])
        for s in current_seq:
            known_f.write(s + "\\n")

with open("${classified_lib}") as fh:
    for line in fh:
        line = line.rstrip("\\n")
        if line.startswith(">"):
            flush()
            current_header = line
            current_seq = []
            # Check classification
            seqid = line.split()[0][1:]  # strip >
            if "#" in seqid:
                classif = seqid.split("#", 1)[1]
                is_unknown = classif.lower().startswith("unknown") or classif == ""
            else:
                is_unknown = True
        else:
            current_seq.append(line)
    flush()

known_f.close()
unknown_f.close()

print(f"Split: {len(known_ids)} known, {len(unknown_ids)} unknown", file=sys.stderr)
PYEOF

    n_unknown=\$(grep -c '^>' unknown_repeats.fa || echo 0)
    n_known=\$(grep -c '^>' known_repeats.fa || echo 0)
    echo "Known families: \${n_known}, Unknown families: \${n_unknown}" >&2

    # ── 2. Run Terrier on Unknown families ────────────────────────────────────
    if [ "\${n_unknown}" -gt 0 ]; then
        terrier \\
            --input unknown_repeats.fa \\
            --output-fasta terrier_classified.fa \\
            --csv terrier_probabilities.csv \\
            --threshold ${threshold}
    else
        touch terrier_classified.fa
        touch terrier_probabilities.csv
        echo "No unknown families to reclassify" >&2
    fi

    # ── 3. Parse Terrier output and build reclassified library ────────────────
    python3 - <<'PYEOF2'
import csv
import sys

THRESHOLD = ${threshold}
ASM = "${assembly_id}"

# ── Parse Terrier classified FASTA ────────────────────────────────────────
#   Terrier appends classification to the header:
#   >famname#Unknown → >famname#Unknown DNA/hAT-Ac
#   or leaves it unchanged if below threshold
terrier_classes = {}  # original_id → new_class

try:
    with open("terrier_classified.fa") as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line.strip()[1:].split()
                original_id = parts[0]  # famname#Unknown
                if len(parts) > 1:
                    # Terrier appended a classification
                    new_class = parts[-1]
                    # Validate it looks like a TE class
                    if "/" in new_class or new_class in (
                        "DNA", "LINE", "SINE", "LTR", "RC",
                        "Retroposon", "Satellite", "Unknown"
                    ):
                        if new_class.lower() != "unknown":
                            terrier_classes[original_id] = new_class
except FileNotFoundError:
    pass

n_reclassified = len(terrier_classes)
print(f"Terrier reclassified {n_reclassified} families (threshold={THRESHOLD})",
      file=sys.stderr)

# ── Build reclassified library ────────────────────────────────────────────
#   Known families: pass through unchanged
#   Unknown families: update header if Terrier provided a classification
n_kept_unknown = 0

with open(f"{ASM}_reclassified_repeats.fa", "w") as out:
    # Write known families unchanged
    with open("known_repeats.fa") as fh:
        for line in fh:
            out.write(line)

    # Write unknown families with updated classification
    current_header = None
    current_seq = []

    def flush_unknown():
        global n_kept_unknown
        if current_header is None:
            return
        seqid = current_header.split()[0][1:]
        rest = current_header.split()[1:] if len(current_header.split()) > 1 else []

        if seqid in terrier_classes:
            # Replace Unknown with Terrier classification
            base_name = seqid.split("#")[0]
            new_class = terrier_classes[seqid]
            new_header = f">{base_name}#{new_class}"
            if rest:
                new_header += " " + " ".join(rest)
            out.write(new_header + "\\n")
        else:
            # Keep as Unknown
            out.write(current_header + "\\n")
            n_kept_unknown += 1

        for s in current_seq:
            out.write(s + "\\n")

    try:
        with open("unknown_repeats.fa") as fh:
            for line in fh:
                line = line.rstrip("\\n")
                if line.startswith(">"):
                    flush_unknown()
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            flush_unknown()
    except FileNotFoundError:
        pass

print(f"Final: {n_reclassified} reclassified, {n_kept_unknown} still Unknown",
      file=sys.stderr)

# ── Write report ──────────────────────────────────────────────────────────
with open(f"{ASM}_reclassification_report.tsv", "w") as out:
    out.write("repeat_id\\toriginal_class\\tterrier_class\\tstatus\\n")

    # Parse all unknown IDs
    try:
        with open("unknown_repeats.fa") as fh:
            for line in fh:
                if line.startswith(">"):
                    seqid = line.strip().split()[0][1:]
                    original = seqid.split("#", 1)[1] if "#" in seqid else "Unknown"
                    if seqid in terrier_classes:
                        out.write(f"{seqid}\\t{original}\\t{terrier_classes[seqid]}\\treclassified\\n")
                    else:
                        out.write(f"{seqid}\\t{original}\\t-\\tunchanged\\n")
    except FileNotFoundError:
        pass

PYEOF2

    # ── 4. Summary ────────────────────────────────────────────────────────────
    n_total=\$(grep -c '^>' ${assembly_id}_reclassified_repeats.fa)
    n_still_unknown=\$(grep '^>' ${assembly_id}_reclassified_repeats.fa | grep -ci '#unknown' || echo 0)
    echo "Reclassified library: \${n_total} total families, \${n_still_unknown} still Unknown" >&2
    """

    stub:
    """
    touch ${assembly_id}_reclassified_repeats.fa
    touch ${assembly_id}_reclassification_report.tsv
    """
}