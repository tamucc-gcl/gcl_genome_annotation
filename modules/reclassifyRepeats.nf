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

    Terrier's --unknown-only flag ensures that already-classified families
    (those with a class/superfamily from RepeatClassifier) pass through
    unchanged — only "Unknown" headers are updated.

    Inputs:
        - assembly_id:       Sample identifier
        - classified_lib:    RepeatClassifier output FASTA (headers have #Class)

    Outputs:
        - reclassified_lib:  ${assembly_id}_reclassified_repeats.fa
        - reclass_report:    ${assembly_id}_reclassification_report.tsv

    Tools: Terrier (pip: bio-terrier), Python 3 (stdlib only)
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
    path terrier_runner

    output:
    tuple val(assembly_id), path("${assembly_id}_reclassified_repeats.fa"),
        emit: reclassified_lib
    tuple val(assembly_id), path("${assembly_id}_reclassification_report.tsv"),
        emit: reclass_report

    script:
    def threshold = params.terrier_threshold == null ? 0.7 : params.terrier_threshold
    """
    set -euo pipefail
    if [ ! -s '${classified_lib}' ]; then
        touch ${assembly_id}_reclassified_repeats.fa
        printf 'repeat_id\\toriginal_class\\tterrier_class\\tstatus\\n' > ${assembly_id}_reclassification_report.tsv
        exit 0
    fi

    # ── 1. Run Terrier on entire library (--unknown-only renames only Unknown) ─
    python3 '${terrier_runner}' \\
        --input ${classified_lib} \\
        --output-fasta ${assembly_id}_reclassified_repeats.fa \\
        --output-csv terrier_probabilities.csv \\
        --threshold ${threshold} \\
        --unknown-only

    # ── 2. Build reclassification report ──────────────────────────────────────
    python3 - <<'PYEOF'
import sys

ASM = "${assembly_id}"

# Compare input and output headers to identify what changed
original_headers = {}
with open("${classified_lib}") as fh:
    for line in fh:
        if line.startswith(">"):
            seqid = line.strip().split()[0][1:]
            original_headers[seqid] = line.strip()

reclassified_headers = {}
with open(f"{ASM}_reclassified_repeats.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            seqid_full = line.strip().split()[0][1:]
            reclassified_headers[seqid_full] = line.strip()

# Match by sequence order (Terrier may change the seqid when reclassifying)
orig_list = []
with open("${classified_lib}") as fh:
    for line in fh:
        if line.startswith(">"):
            orig_list.append(line.strip())

reclass_list = []
with open(f"{ASM}_reclassified_repeats.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            reclass_list.append(line.strip())

n_reclassified = 0
n_unchanged = 0
n_known = 0
if len(orig_list) != len(reclass_list):
    raise ValueError('Terrier changed the number of repeat families')

with open(f"{ASM}_reclassification_report.tsv", "w") as out:
    out.write("repeat_id\\toriginal_class\\tterrier_class\\tstatus\\n")

    for i, (orig_hdr, new_hdr) in enumerate(zip(orig_list, reclass_list)):
        orig_id = orig_hdr.split()[0][1:]
        new_id = new_hdr.split()[0][1:]

        orig_class = orig_id.split("#", 1)[1] if "#" in orig_id else ""
        new_class = new_id.split("#", 1)[1] if "#" in new_id else ""

        is_unknown = orig_class.lower().startswith("unknown") or orig_class == ""

        if not is_unknown:
            n_known += 1
            # Don't report known families
            continue

        if new_class != orig_class and new_class.lower() != "unknown" and new_class != "":
            n_reclassified += 1
            out.write(f"{orig_id}\\t{orig_class}\\t{new_class}\\treclassified\\n")
        else:
            n_unchanged += 1
            out.write(f"{orig_id}\\t{orig_class}\\t-\\tunchanged\\n")

print(f"Reclassification: {n_known} already classified, "
      f"{n_reclassified} reclassified by Terrier, "
      f"{n_unchanged} still Unknown", file=sys.stderr)
PYEOF

    # ── 3. Summary ────────────────────────────────────────────────────────────
    n_total=\$(grep -c '^>' ${assembly_id}_reclassified_repeats.fa)
    n_still_unknown=\$(grep '^>' ${assembly_id}_reclassified_repeats.fa | grep -ci '#unknown' || true)
    echo "Reclassified library: \${n_total} total families, \${n_still_unknown} still Unknown" >&2
    """

    stub:
    """
    touch ${assembly_id}_reclassified_repeats.fa
    touch ${assembly_id}_reclassification_report.tsv
    """
}
