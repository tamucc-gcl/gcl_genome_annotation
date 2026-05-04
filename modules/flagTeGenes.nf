/*
========================================================================================
    FLAG TE-DERIVED GENES
========================================================================================
    Identify and flag predicted gene models that are likely derived from
    transposable elements rather than bona fide host genes.

    Evidence sources:
      1. Genomic overlap: gene model is ≥50% covered by RepeatMasker annotations
      2. Protein homology: predicted protein has significant BLASTP similarity
         to known TE proteins (from the de novo repeat library, translated)

    Flagged genes receive a Note=putative_TE_derived attribute in the GFF3.
    They are NOT removed — downstream users can filter at their discretion.

    Also produces a summary TSV with per-gene evidence for the flag.

    Inputs:
        - assembly_id:       Sample identifier
        - braker_gff3:       BRAKER3 gene models GFF3
        - braker_proteins:   BRAKER3 predicted proteins FASTA
        - repeat_gff:        RepeatMasker GFF annotation
        - repeat_lib:        De novo repeat library FASTA (pre-filter, includes
                             all TE sequences for homology search)

    Outputs:
        - flagged_gff3:      ${assembly_id}.te_flagged.gff3
        - te_gene_report:    ${assembly_id}.te_gene_report.tsv

    Tools: bedtools, DIAMOND (blastp, makedb), Python 3 (stdlib only)
========================================================================================
*/

process FLAG_TE_GENES {
    tag "${assembly_id}"
    label 'flag_te_genes'
    publishDir "${params.outdir}/structural_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(braker_gff3), path(braker_proteins),
        path(repeat_gff), path(repeat_lib)

    output:
    tuple val(assembly_id), path("${assembly_id}.te_flagged.gff3"),
        emit: flagged_gff3
    tuple val(assembly_id), path("${assembly_id}.te_gene_report.tsv"),
        emit: te_gene_report

    script:
    def overlap_frac = params.te_flag_overlap_frac ?: 0.5
    def blastp_evalue = params.te_flag_evalue ?: 1e-10
    """
    set -euo pipefail

    # ── 1. Find genes with ≥N% overlap with repeat annotations ────────────────
    #   Extract gene features, sort, intersect with repeat GFF
    awk '\$3 == "gene"' ${braker_gff3} | sort -k1,1 -k4,4n > genes.bed

    # bedtools intersect: -f = min fraction of gene covered by repeats
    # -wo = write overlap bp
    bedtools intersect \
        -a genes.bed \
        -b <(awk '!/^#/ && NF >= 9' ${repeat_gff} | sort -k1,1 -k4,4n) \
        -f ${overlap_frac} \
        -wo \
        | awk '{
            # Extract gene ID from attributes (col 9)
            match(\$9, /ID=([^;]+)/, arr)
            if (arr[1] != "") print arr[1]
        }' \
        | sort -u > overlap_flagged_ids.txt

    echo "Genes with >=${overlap_frac} repeat overlap: \$(wc -l < overlap_flagged_ids.txt)" >&2

    # ── 2. Translate repeat library to protein and build DIAMOND DB ───────────
    #   Use all 6 frames to capture TE ORFs
    seqkit translate -f 6 --min-len 30 ${repeat_lib} > repeat_lib_translated.fa

    diamond makedb \
        --in repeat_lib_translated.fa \
        --db repeat_lib_db \
        --threads ${task.cpus}

    # ── 3. BLASTP predicted proteins against translated repeat library ────────
    diamond blastp \
        --db repeat_lib_db \
        --query ${braker_proteins} \
        --out prot_vs_repeats.tsv \
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
        --evalue ${blastp_evalue} \
        --max-target-seqs 1 \
        --threads ${task.cpus} \
        || true  # OK if no hits

    # ── 4. Identify proteins with strong TE homology ──────────────────────────
    #   Map protein IDs back to gene IDs and flag
    python3 - <<'PYEOF'
import re
import sys
from collections import defaultdict

ASM = "${assembly_id}"
OVERLAP_FRAC = ${overlap_frac}

# ── Load overlap-flagged gene IDs ─────────────────────────────────────────
overlap_genes = set()
with open("overlap_flagged_ids.txt") as fh:
    for line in fh:
        overlap_genes.add(line.strip())

# ── Parse BRAKER GFF3 to build mRNA→gene mapping ─────────────────────────
mrna_to_gene = {}
gene_ids = set()

with open("${braker_gff3}") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.strip().split("\\t")
        if len(cols) < 9:
            continue
        ftype = cols[2]
        attrs = {}
        for token in cols[8].split(";"):
            if "=" in token:
                k, v = token.split("=", 1)
                attrs[k.strip()] = v.strip()

        if ftype == "gene":
            gene_ids.add(attrs.get("ID", ""))
        elif ftype in ("mRNA", "transcript"):
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            mrna_to_gene[fid] = parent

# ── Parse BLASTP results ──────────────────────────────────────────────────
#   Protein IDs from BRAKER are typically mRNA IDs (gene_id.t1)
blastp_genes = {}  # gene_id → best hit info

with open("prot_vs_repeats.tsv") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.strip().split("\\t")
        if len(cols) < 8:
            continue
        prot_id = cols[0]
        hit_id  = cols[1]
        pident  = float(cols[2])
        alnlen  = int(cols[3])
        qlen    = int(cols[4])
        evalue  = float(cols[6])
        bitscore = float(cols[7])

        # Query coverage
        qcov = (alnlen / qlen * 100) if qlen > 0 else 0

        # Map protein ID to gene ID
        gene_id = mrna_to_gene.get(prot_id, "")
        if not gene_id:
            # Try stripping .t1 suffix
            base = prot_id.rsplit(".", 1)[0] if "." in prot_id else prot_id
            if base in gene_ids:
                gene_id = base
            else:
                gene_id = mrna_to_gene.get(base, "")

        if not gene_id:
            continue

        # Keep best hit per gene
        if gene_id not in blastp_genes or bitscore > blastp_genes[gene_id]["bitscore"]:
            blastp_genes[gene_id] = {
                "prot_id": prot_id, "hit_id": hit_id,
                "pident": pident, "qcov": qcov,
                "evalue": evalue, "bitscore": bitscore,
            }

blastp_flagged = set(blastp_genes.keys())

# ── Combine evidence ──────────────────────────────────────────────────────
all_flagged = overlap_genes | blastp_flagged

# Evidence categories
both = overlap_genes & blastp_flagged
overlap_only = overlap_genes - blastp_flagged
blastp_only = blastp_flagged - overlap_genes

print(f"TE-derived gene flags:", file=sys.stderr)
print(f"  Repeat overlap only: {len(overlap_only)}", file=sys.stderr)
print(f"  Protein homology only: {len(blastp_only)}", file=sys.stderr)
print(f"  Both: {len(both)}", file=sys.stderr)
print(f"  Total flagged: {len(all_flagged)}", file=sys.stderr)

# ── Write gene-level report ───────────────────────────────────────────────
with open(f"{ASM}.te_gene_report.tsv", "w") as out:
    out.write("gene_id\\tflagged\\tevidence\\trepeat_overlap\\t"
              "blastp_hit\\tbp_pident\\tbp_qcov\\tbp_evalue\\n")
    for gid in sorted(gene_ids):
        flagged = gid in all_flagged
        evidence = []
        if gid in overlap_genes:
            evidence.append("repeat_overlap")
        if gid in blastp_flagged:
            evidence.append("protein_homology")

        has_overlap = "yes" if gid in overlap_genes else "no"

        bp = blastp_genes.get(gid, {})
        bp_hit = bp.get("hit_id", "-")
        bp_pident = f"{bp['pident']:.1f}" if bp else "-"
        bp_qcov = f"{bp['qcov']:.1f}" if bp else "-"
        bp_evalue = f"{bp['evalue']:.2e}" if bp else "-"

        out.write(f"{gid}\\t{'yes' if flagged else 'no'}\\t"
                  f"{';'.join(evidence) if evidence else '-'}\\t"
                  f"{has_overlap}\\t{bp_hit}\\t{bp_pident}\\t{bp_qcov}\\t{bp_evalue}\\n")

# ── Write flagged gene IDs for GFF3 annotation ───────────────────────────
with open("all_flagged_ids.txt", "w") as out:
    for gid in sorted(all_flagged):
        out.write(gid + "\\n")

# ── Build evidence strings per gene for GFF3 Note ─────────────────────────
with open("flag_notes.tsv", "w") as out:
    for gid in sorted(all_flagged):
        parts = []
        if gid in overlap_genes:
            parts.append(f"repeat_overlap>={OVERLAP_FRAC:.0%}")
        if gid in blastp_genes:
            bp = blastp_genes[gid]
            parts.append(f"TE_homology:{bp['hit_id']}:{bp['pident']:.0f}%id")
        out.write(f"{gid}\\t{';'.join(parts)}\\n")

PYEOF

    # ── 5. Annotate GFF3 with TE flags ────────────────────────────────────────
    python3 - <<'PYEOF2'
import sys
from collections import OrderedDict

# Load flag notes
flag_notes = {}
with open("flag_notes.tsv") as fh:
    for line in fh:
        parts = line.strip().split("\\t", 1)
        if len(parts) == 2:
            flag_notes[parts[0]] = parts[1]

# Load all flagged IDs (genes)
flagged_ids = set()
with open("all_flagged_ids.txt") as fh:
    for line in fh:
        flagged_ids.add(line.strip())

# Build gene→flag mapping, then propagate to children
# We need to know which genes are flagged, then flag all their children too
gene_to_children = {}  # gene_id → set of descendant IDs
parent_map = {}        # feature_id → parent_id

# First pass: build parent map
features = []
with open("${braker_gff3}") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            features.append(("comment", line))
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) != 9:
            features.append(("comment", line))
            continue

        attrs = OrderedDict()
        for token in cols[8].split(";"):
            token = token.strip()
            if "=" in token:
                k, v = token.split("=", 1)
                attrs[k] = v

        fid = attrs.get("ID", "")
        parent = attrs.get("Parent", "")
        if parent:
            parent_map[fid] = parent

        features.append(("feature", cols[:8], attrs))

# Propagate: find all features whose ancestor gene is flagged
def get_gene_ancestor(fid):
    \"\"\"Walk up Parent chain to find gene-level ancestor.\"\"\"
    visited = set()
    current = fid
    while current in parent_map and current not in visited:
        visited.add(current)
        current = parent_map[current]
    return current

flagged_features = set()
for tag, *rest in features:
    if tag != "feature":
        continue
    cols8, attrs = rest
    fid = attrs.get("ID", "")
    ancestor = get_gene_ancestor(fid)
    if ancestor in flagged_ids:
        flagged_features.add(fid)
        flagged_features.add(ancestor)

# Second pass: write output with flags
with open("${assembly_id}.te_flagged.gff3", "w") as out:
    for item in features:
        if item[0] == "comment":
            out.write(item[1])
            continue

        cols8, attrs = item[1], item[2]
        fid = attrs.get("ID", "")

        if fid in flagged_features:
            # Add/append TE flag to gene-level features
            ancestor = get_gene_ancestor(fid)
            note_text = flag_notes.get(ancestor, "TE_evidence")

            if cols8[2] == "gene":
                existing = attrs.get("Note", "")
                new_note = f"putative_TE_derived({note_text})"
                attrs["Note"] = f"{existing};{new_note}" if existing else new_note

            # Add simpler flag to all child features
            elif fid in flagged_features:
                existing = attrs.get("Note", "")
                flag = "parent_flagged_TE_derived"
                attrs["Note"] = f"{existing};{flag}" if existing else flag

        # Serialize attributes
        attr_str = ";".join(
            f"{k}={v}" if v else k for k, v in attrs.items()
        )
        out.write("\\t".join(cols8 + [attr_str]) + "\\n")

n_flagged_genes = len(flagged_ids)
n_total_genes = sum(1 for item in features
                    if item[0] == "feature" and item[1][2] == "gene")
print(f"Flagged {n_flagged_genes}/{n_total_genes} genes as putative TE-derived",
      file=sys.stderr)

PYEOF2
    """

    stub:
    """
    touch ${assembly_id}.te_flagged.gff3
    touch ${assembly_id}.te_gene_report.tsv
    """
}