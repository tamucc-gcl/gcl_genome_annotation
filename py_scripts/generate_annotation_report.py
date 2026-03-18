#!/usr/bin/env python3
"""
generate_annotation_report.py
=============================
Unified multi-assembly Markdown annotation report.
Cross-assembly comparison tables with per-assembly detail where needed.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime


def img(pd, filename):
    path = os.path.join(pd, filename)
    if os.path.exists(path):
        return f"![{filename}]({pd}/{filename})\n"
    return ""


def fmt(n):
    try:
        return f"{int(float(n)):,}"
    except (ValueError, TypeError):
        return str(n)


def read_text(path):
    try:
        with open(path) as f:
            return f.read()
    except (FileNotFoundError, IOError):
        return ""


def infer_id(path, suffixes):
    base = os.path.basename(path)
    for s in suffixes:
        if base.endswith(s):
            return base[: -len(s)]
    return os.path.splitext(base)[0]


def find_file(file_list, assembly_id, suffixes):
    for f in file_list:
        if infer_id(f, suffixes) == assembly_id:
            return f
    return None


def extract_omark_pct(text, key):
    m = re.search(re.escape(key) + r"([\d.]+)%", text)
    return m.group(1) if m else "—"


def parse_omark_summary(text):
    result = {}
    for line in text.splitlines():
        ls = line.strip().lstrip("#").strip()
        if ls.startswith("S:") and "%" in ls:
            result["comp_line"] = ls
        if ls.startswith("A:") and "%" in ls:
            result["cons_line"] = ls
        if "selected clade was" in ls:
            result["clade"] = ls.split("was")[-1].strip()
        if "conserved HOGs" in ls.lower():
            m = re.search(r"(\d[\d,]*)", ls)
            if m:
                result["n_hogs"] = m.group(1)
        if "there are" in ls.lower() and "proteins" in ls.lower():
            m = re.search(r"(\d[\d,]*)\s+proteins", ls)
            if m:
                result["n_proteins"] = m.group(1)
    return result


def multi_col_header(assembly_ids):
    """Markdown table header with assembly columns."""
    h = "| Metric | " + " | ".join(f"`{a}`" for a in assembly_ids) + " |"
    s = "|:---|" + "|".join("---:" for _ in assembly_ids) + "|"
    return h + "\n" + s + "\n"


def multi_col_row(label, rows, key):
    vals = " | ".join(rows[a].get(key, "—") for a in sorted(rows))
    return f"| {label} | {vals} |\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--assembly_stats", nargs="+", required=True)
    parser.add_argument("--repeat_tbls", nargs="+", required=True)
    parser.add_argument("--busco_txts", nargs="+", required=True)
    parser.add_argument("--omark_summaries", nargs="+", required=True)
    parser.add_argument("--omark_detaileds", nargs="+", required=True)
    parser.add_argument("--agat_stats", nargs="+", required=True)
    parser.add_argument("--plots_dir", required=True)
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    pd = args.plots_dir

    # Read metrics keyed by assembly_id
    met = {}
    with open(args.metrics) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            met[row["assembly_id"]] = row
    aids = sorted(met.keys())
    n_asm = len(aids)
    multi = n_asm > 1

    md = []
    md.append("# Genome Annotation Report\n")
    md.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in aids)}\n")
    md.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    md.append("---\n")

    # ── TOC ───────────────────────────────────────────────────────────────
    md.append("""## Table of Contents

1. [Overview](#1-overview)
2. [Gene Model Structure](#2-gene-model-structure)
3. [Repeat Annotation](#3-repeat-annotation)
4. [Non-coding RNA](#4-non-coding-rna)
5. [Telomeres & Centromeres](#5-telomeres--centromeres)
6. [Quality Assessment](#6-quality-assessment)
7. [Functional Annotation](#7-functional-annotation)
8. [Tuning Recommendations](#8-tuning-recommendations)
""")

    # ═════════════════════════════════════════════════════════════════════
    # 1. Overview
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 1. Overview\n")

    md.append(multi_col_header(aids))
    overview = [
        ("n_genes", "Protein-coding genes"),
        ("n_mrnas", "mRNA transcripts"),
        ("n_proteins", "Predicted proteins"),
        ("n_trna", "tRNA genes"),
        ("n_ncrna", "ncRNA features"),
        ("n_telomere", "Telomeric regions"),
        ("n_centromere", "Centromere features"),
        ("n_repeat", "Repeat features"),
        ("scaffolds_with_genes", "Scaffolds with genes"),
        ("repeat_masked_pct", "Genome masked (%)"),
        ("busco_complete", "BUSCO complete (%)"),
        ("omark_consistent", "OMArk consistent (%)"),
        ("n_named", "Named products"),
        ("n_hypothetical", "Hypothetical proteins"),
        ("pct_named", "Named (%)"),
        ("unique_go_terms", "Unique GO terms"),
        ("unique_kegg_kos", "Unique KEGG KOs"),
    ]
    for key, label in overview:
        is_pct = "pct" in key or "busco" in key or "omark" in key or "masked_pct" in key
        vals = " | ".join(
            met[a].get(key, "—") if is_pct else fmt(met[a].get(key, "—"))
            for a in aids
        )
        md.append(f"| {label} | {vals} |\n")
    md.append("\n")

    if multi:
        md.append(img(pd, "C_feature_counts.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 2. Gene Model Structure
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 2. Gene Model Structure\n")

    for aid in aids:
        md.append(f"### {aid}\n") if multi else None
        md.append(img(pd, f"{aid}_gene_length_dist.png"))
        md.append(img(pd, f"{aid}_exons_per_tx.png"))
        md.append(img(pd, f"{aid}_intron_qc.png"))
        md.append(img(pd, f"{aid}_genes_per_scaffold.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 3. Repeat Annotation
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 3. Repeat Annotation\n")

    # Summary table across assemblies
    md.append(multi_col_header(aids))
    md.append(
        "| Genome masked (%) | " +
        " | ".join(met[a].get("repeat_masked_pct", "—") for a in aids) + " |\n"
    )
    md.append(
        "| Masked bp | " +
        " | ".join(fmt(met[a].get("repeat_masked_bp", "—")) for a in aids) + " |\n"
    )
    md.append(
        "| Repeat features | " +
        " | ".join(fmt(met[a].get("n_repeat", "—")) for a in aids) + " |\n"
    )
    md.append("\n")

    md.append(img(pd, "C_repeat_masking.png"))

    for aid in aids:
        tbl_file = find_file(args.repeat_tbls, aid, [".tbl"])
        if not tbl_file:
            continue

        md.append(f"### {aid} — Repeat summary\n") if multi else md.append("### Summary\n")

        tbl_text = read_text(tbl_file).strip()
        # Extract the compact header (through "bases masked" line)
        summary_lines = []
        for line in tbl_text.splitlines():
            summary_lines.append(line)
            if "bases masked" in line.lower():
                break
        if summary_lines:
            md.append("```\n" + "\n".join(summary_lines) + "\n```\n")

        # Per-scaffold masking percentage
        # Parse .tbl for per-sequence stats if available, otherwise
        # note that the full report is available.
        md.append(
            f"> Full repeat tables: "
            f"`results/repeats/{aid}.tbl`, "
            f"`results/repeats/{aid}.buildSummary.tbl`, "
            f"`results/repeats/{aid}.repeat_landscape.html`\n\n"
        )

    # ═════════════════════════════════════════════════════════════════════
    # 4. Non-coding RNA
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 4. Non-coding RNA\n")

    md.append(multi_col_header(aids))
    md.append("| tRNA genes | " + " | ".join(fmt(met[a].get("n_trna", "—")) for a in aids) + " |\n")
    md.append("| ncRNA features (Rfam) | " + " | ".join(fmt(met[a].get("n_ncrna", "—")) for a in aids) + " |\n")
    md.append("\n")

    for aid in aids:
        md.append(f"### {aid}\n") if multi else None
        md.append(img(pd, f"{aid}_trna_isotypes.png"))
        md.append(img(pd, f"{aid}_ncrna_rfam.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 5. Telomeres & Centromeres
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 5. Telomeres & Centromeres\n")

    md.append(multi_col_header(aids))
    md.append("| Telomeric regions | " + " | ".join(fmt(met[a].get("n_telomere", "—")) for a in aids) + " |\n")
    md.append("| Centromere features | " + " | ".join(fmt(met[a].get("n_centromere", "—")) for a in aids) + " |\n")
    md.append("\n")

    for aid in aids:
        md.append(f"### {aid}\n") if multi else None
        md.append(img(pd, f"{aid}_telomere_positions.png"))
        md.append(img(pd, f"{aid}_telomere_density.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 6. Quality Assessment
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 6. Quality Assessment\n")

    # ── BUSCO table ───────────────────────────────────────────────────
    md.append("### BUSCO\n")
    md.append(multi_col_header(aids))
    busco_fields = [
        ("busco_complete", "Complete (%)"),
        ("busco_single", "— Single (%)"),
        ("busco_duplicated", "— Duplicated (%)"),
        ("busco_fragmented", "Fragmented (%)"),
        ("busco_missing", "Missing (%)"),
        ("busco_total", "Total BUSCOs"),
    ]
    for key, label in busco_fields:
        vals = " | ".join(met[a].get(key, "—") for a in aids)
        md.append(f"| {label} | {vals} |\n")

    # Extract lineage from first BUSCO file
    for aid in aids:
        bf = find_file(args.busco_txts, aid, ["_busco.txt"])
        if bf:
            for line in read_text(bf).splitlines():
                if "lineage dataset" in line.lower():
                    md.append(f"\n*{line.strip()}*\n")
                    break
            break
    md.append("\n")

    md.append(img(pd, "C_busco.png"))

    md.append(f"> Full BUSCO reports: `results/qc/busco/`\n\n")

    # ── OMArk tables ──────────────────────────────────────────────────
    md.append("### OMArk\n")

    # Parse per-assembly OMArk summaries for clade/HOG/protein counts
    omark_info = {}
    for aid in aids:
        of = find_file(args.omark_summaries, aid, [".summary.txt"])
        if of:
            omark_info[aid] = parse_omark_summary(read_text(of))

    # Show clade and HOG count
    if omark_info:
        first = next(iter(omark_info.values()))
        md.append(f"**Reference clade:** {first.get('clade', '—')}\n\n")

    md.append("#### Completeness (conserved genes)\n")
    md.append(multi_col_header(aids))
    md.append("| Conserved HOGs | " + " | ".join(
        omark_info.get(a, {}).get("n_hogs", "—") for a in aids) + " |\n")
    omark_comp = [
        ("omark_single", "Single (%)"),
        ("omark_duplicated", "Duplicated total (%)"),
        ("omark_dup_unexpected", "— Unexpected (%)"),
        ("omark_dup_expected", "— Expected (%)"),
        ("omark_missing", "Missing (%)"),
    ]
    for key, label in omark_comp:
        vals = " | ".join(met[a].get(key, "—") for a in aids)
        md.append(f"| {label} | {vals} |\n")
    md.append("\n")

    md.append(img(pd, "C_omark_completeness.png"))

    md.append("#### Consistency (proteome quality)\n")
    md.append(multi_col_header(aids))
    md.append("| Total proteins | " + " | ".join(
        omark_info.get(a, {}).get("n_proteins", "—") for a in aids) + " |\n")
    omark_cons = [
        ("omark_consistent", "Consistent (%)"),
        ("omark_partial", "— Partial hits (%)"),
        ("omark_fragmented", "— Fragmented (%)"),
        ("omark_inconsistent", "Inconsistent (%)"),
        ("omark_contaminant", "Contaminant (%)"),
        ("omark_unknown", "Unknown (%)"),
    ]
    for key, label in omark_cons:
        vals = " | ".join(met[a].get(key, "—") for a in aids)
        md.append(f"| {label} | {vals} |\n")
    md.append("\n")

    md.append(img(pd, "C_omark_consistency.png"))

    md.append(f"> Full OMArk reports: `results/qc/omark/`\n\n")

    # ── AGAT summary table ────────────────────────────────────────────
    md.append("### Structural statistics (AGAT)\n")

    # Parse AGAT stats for each assembly into a dict
    agat_data = {}
    agat_keys_ordered = []
    key_patterns = {
        "number of gene": "Genes",
        "number of mrna": "mRNAs",
        "number of exon": "Exons",
        "number of cds": "CDS features",
        "number of intron": "Introns",
        "number of single exon gene": "Single-exon genes",
        "number of single exon mrna": "Single-exon mRNAs",
        "mean gene length": "Mean gene length (bp)",
        "mean mrna length": "Mean mRNA length (bp)",
        "mean cds length": "Mean CDS length (bp)",
        "mean exon length": "Mean exon length (bp)",
        "mean intron length": "Mean intron length (bp)",
        "longest gene": "Longest gene (bp)",
        "shortest gene": "Shortest gene (bp)",
        "total gene length": "Total gene space (bp)",
        "mean exons per": "Mean exons per mRNA",
        "mean introns per": "Mean introns per mRNA",
        "% of genome": "% genome covered by genes",
    }

    for aid in aids:
        agat_file = find_file(args.agat_stats, aid, ["_annotation_stats"])
        if not agat_file:
            continue
        agat_text = read_text(agat_file)
        agat_data[aid] = {}
        for line in agat_text.splitlines():
            ll = line.lower().strip()
            for pattern, label in key_patterns.items():
                if pattern in ll:
                    # Extract the numeric value (last number-like token on line)
                    nums = re.findall(r"[\d,.]+", line)
                    if nums:
                        val = nums[-1]
                        agat_data[aid][label] = val
                        if label not in agat_keys_ordered:
                            agat_keys_ordered.append(label)

    if agat_data:
        md.append(multi_col_header(aids))
        for label in agat_keys_ordered:
            vals = " | ".join(agat_data.get(a, {}).get(label, "—") for a in aids)
            md.append(f"| {label} | {vals} |\n")
        md.append("\n")

    md.append(f"> Full AGAT reports: `results/qc/`\n\n")

    if multi:
        md.append(img(pd, "C_quality_summary.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 7. Functional Annotation
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 7. Functional Annotation\n")

    md.append("### Tool coverage\n")
    md.append(multi_col_header(aids))
    func_fields = [
        ("n_proteins", "Total proteins"),
        ("n_ipr", "InterProScan hits"),
        ("pct_ipr", "InterProScan (%)"),
        ("n_eggnog", "eggNOG hits"),
        ("pct_eggnog", "eggNOG (%)"),
        ("n_diamond", "DIAMOND hits"),
        ("pct_diamond", "DIAMOND (%)"),
        ("n_kofamscan", "KofamScan hits"),
        ("pct_kofamscan", "KofamScan (%)"),
        ("n_deepsig", "DeepSig hits"),
        ("pct_deepsig", "DeepSig (%)"),
        ("n_named", "Named products"),
        ("n_hypothetical", "Hypothetical"),
        ("pct_named", "Named (%)"),
        ("pct_hypothetical", "Hypothetical (%)"),
        ("unique_go_terms", "Unique GO terms"),
        ("unique_kegg_kos", "Unique KEGG KOs"),
    ]
    for key, label in func_fields:
        is_pct = "pct" in key
        vals = " | ".join(
            met[a].get(key, "—") if is_pct else fmt(met[a].get(key, "—"))
            for a in aids
        )
        md.append(f"| {label} | {vals} |\n")
    md.append("\n")

    if multi:
        md.append(img(pd, "C_functional_coverage.png"))
        md.append(img(pd, "C_product_names.png"))
        md.append(img(pd, "C_go_kegg_diversity.png"))

    for aid in aids:
        md.append(f"### {aid}\n") if multi else md.append("### Details\n")
        md.append(img(pd, f"{aid}_functional_coverage.png"))
        md.append(img(pd, f"{aid}_top_products.png"))
        md.append(img(pd, f"{aid}_top_go.png"))
        md.append(img(pd, f"{aid}_go_treemap.png"))
        md.append(img(pd, f"{aid}_go_ontology_composition.png"))
        md.append(img(pd, f"{aid}_top_kegg.png"))
        md.append(img(pd, f"{aid}_top_ipr.png"))

    # ═════════════════════════════════════════════════════════════════════
    # 8. Tuning Recommendations
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 8. Tuning Recommendations\n")

    any_flags = False
    for aid in aids:
        row = met[aid]
        flags = []

        try:
            v = float(row.get("pct_hypothetical", 0))
            if v > 50:
                flags.append(f"High hypothetical rate ({v:.1f}%) — consider RNA-seq "
                             "evidence for BRAKER3 or broader protein databases.")
        except ValueError:
            pass

        try:
            v = float(row.get("busco_complete", 0))
            if v < 80:
                flags.append(f"Low BUSCO completeness ({v:.1f}%) — check repeat masking, "
                             "assembly contiguity, or increase BRAKER3 `--rounds`.")
        except ValueError:
            pass

        try:
            v = float(row.get("omark_unknown", 0))
            if v > 25:
                flags.append(f"High OMArk unknown rate ({v:.1f}%) — possible over-prediction. "
                             "Consider TSEBRA filtering.")
        except ValueError:
            pass

        try:
            v = float(row.get("omark_fragmented", 0))
            if v > 10:
                flags.append(f"High OMArk fragmented rate ({v:.1f}%) — many genes shorter "
                             "than expected. Check assembly gaps.")
        except ValueError:
            pass

        try:
            v = float(row.get("omark_dup_unexpected", 0))
            if v > 15:
                flags.append(f"High unexpected duplication ({v:.1f}%) — possible residual "
                             "haplotig duplication. Cross-reference BUSCO D% and Hi-C maps.")
        except ValueError:
            pass

        try:
            v = int(float(row.get("short_introns_lt10", 0)))
            if v > 50:
                flags.append(f"{v:,} introns < 10 bp — likely gene model artifacts. "
                             "Consider AGAT/TSEBRA post-filtering.")
        except ValueError:
            pass

        try:
            v = float(row.get("repeat_masked_pct", 0))
            if v < 20:
                flags.append(f"Low repeat masking ({v:.1f}%) — fish genomes typically "
                             "30-50%. Try more RepeatModeler rounds or slower RM speed.")
        except ValueError:
            pass

        try:
            v = float(row.get("pct_diamond", 0))
            if v < 30:
                flags.append(f"Low SwissProt coverage ({v:.1f}%) — consider TrEMBL/NR "
                             "for broader homology.")
        except ValueError:
            pass

        if flags:
            any_flags = True
            md.append(f"### `{aid}`\n")
            for f in flags:
                md.append(f"- {f}\n")
            md.append("\n")

    if not any_flags:
        md.append("No major issues detected. Annotation parameters appear well-tuned.\n")

    md.append("---\n*Report generated by the genome annotation pipeline.*\n")

    with open(args.outfile, "w") as out:
        out.write("\n".join(md))

    print(f"Report: {n_asm} assemblies → {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()