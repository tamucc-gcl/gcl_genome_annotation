#!/usr/bin/env python3
"""
generate_annotation_report.py
=============================
Assemble a unified multi-assembly Markdown annotation report.
Includes per-assembly detail sections and cross-assembly comparison.
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter
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


def pct_str(n, total):
    try:
        n, total = float(n), float(total)
    except (ValueError, TypeError):
        return "—"
    if total == 0:
        return "0.0%"
    return f"{100 * n / total:.1f}%"


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
    """Find the file for a given assembly_id from a list of paths."""
    for f in file_list:
        if infer_id(f, suffixes) == assembly_id:
            return f
    return None


def parse_omark_summary(text):
    """Parse OMArk compact summary format."""
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


def extract_omark_pct(text, key):
    m = re.search(re.escape(key) + r"([\d.]+)%", text)
    return m.group(1) if m else "—"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--assembly_stats", nargs="+", required=True)
    parser.add_argument("--repeat_tbls", nargs="+", required=True)
    parser.add_argument("--repeat_tbls_full", nargs="+", required=True)
    parser.add_argument("--busco_txts", nargs="+", required=True)
    parser.add_argument("--omark_summaries", nargs="+", required=True)
    parser.add_argument("--omark_detaileds", nargs="+", required=True)
    parser.add_argument("--agat_stats", nargs="+", required=True)
    parser.add_argument("--plots_dir", required=True)
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    pd = args.plots_dir

    # Read metrics
    rows = []
    with open(args.metrics) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    assembly_ids = [r["assembly_id"] for r in rows]
    n_asm = len(assembly_ids)
    multi = n_asm > 1

    md = []
    md.append("# Genome Annotation Report\n")
    md.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in assembly_ids)}\n")
    md.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    md.append("---\n")

    # ── TOC ───────────────────────────────────────────────────────────────
    md.append("## Table of Contents\n")
    toc = [
        "1. [Overview](#1-overview)",
    ]
    if multi:
        toc.append("2. [Cross-Assembly Comparison](#2-cross-assembly-comparison)")
        toc.append("3. [Per-Assembly Details](#3-per-assembly-details)")
        toc.append("4. [Tuning Recommendations](#4-tuning-recommendations)")
    else:
        toc.append("2. [Assembly Detail](#2-assembly-detail)")
        toc.append("3. [Tuning Recommendations](#3-tuning-recommendations)")
    md.append("\n".join(toc) + "\n")

    # ═════════════════════════════════════════════════════════════════════
    # 1. Overview
    # ═════════════════════════════════════════════════════════════════════
    md.append("---\n## 1. Overview\n")

    overview_fields = [
        ("n_genes", "Protein-coding genes"),
        ("n_mrnas", "mRNA transcripts"),
        ("n_trna", "tRNA genes"),
        ("n_ncrna", "ncRNA features"),
        ("n_telomere", "Telomeric regions"),
        ("n_centromere", "Centromere features"),
        ("n_repeat", "Repeat features"),
        ("n_proteins", "Predicted proteins"),
        ("n_named", "Named products"),
        ("n_hypothetical", "Hypothetical proteins"),
        ("pct_named", "Named (%)"),
        ("repeat_masked_pct", "Repeat masked (%)"),
        ("busco_complete", "BUSCO complete (%)"),
        ("omark_consistent", "OMArk consistent (%)"),
        ("unique_go_terms", "Unique GO terms"),
        ("unique_kegg_kos", "Unique KEGG KOs"),
    ]

    header = "| Metric | " + " | ".join(f"`{a}`" for a in assembly_ids) + " |"
    sep = "|:---|" + "|".join("---:" for _ in assembly_ids) + "|"
    md.append(header + "\n" + sep + "\n")
    for key, label in overview_fields:
        vals = " | ".join(fmt(r.get(key, "—")) if "pct" not in key and "busco" not in key and "omark" not in key
                          else r.get(key, "—")
                          for r in rows)
        md.append(f"| {label} | {vals} |\n")
    md.append("\n")

    # ═════════════════════════════════════════════════════════════════════
    # 2. Cross-Assembly Comparison (only if multi)
    # ═════════════════════════════════════════════════════════════════════
    if multi:
        md.append("---\n## 2. Cross-Assembly Comparison\n")

        md.append("### Feature counts\n")
        md.append(img(pd, "C_feature_counts.png"))

        md.append("### BUSCO completeness\n")
        md.append(img(pd, "C_busco.png"))

        md.append("### OMArk — conserved gene completeness\n")
        md.append(img(pd, "C_omark_completeness.png"))

        md.append("### OMArk — proteome consistency\n")
        md.append(img(pd, "C_omark_consistency.png"))

        md.append("### Functional tool coverage\n")
        md.append(img(pd, "C_functional_coverage.png"))

        md.append("### Product name assignment\n")
        md.append(img(pd, "C_product_names.png"))

        md.append("### Functional vocabulary (GO / KEGG)\n")
        md.append(img(pd, "C_go_kegg_diversity.png"))

        md.append("### Repeat masking\n")
        md.append(img(pd, "C_repeat_masking.png"))

        md.append("### Quality summary\n")
        md.append(img(pd, "C_quality_summary.png"))

    # ═════════════════════════════════════════════════════════════════════
    # Per-Assembly Details
    # ═════════════════════════════════════════════════════════════════════
    section_num = 3 if multi else 2
    md.append(f"---\n## {section_num}. {'Per-Assembly Details' if multi else 'Assembly Detail'}\n")

    for row in rows:
        aid = row["assembly_id"]
        md.append(f"---\n### {aid}\n")

        # ── Assembly stats ────────────────────────────────────────────
        stats_file = find_file(args.assembly_stats, aid,
                               [".assembly_stats.tsv"])
        if stats_file:
            md.append("#### Assembly statistics\n")
            md.append("```\n" + read_text(stats_file).strip() + "\n```\n")

        # ── Gene model structure ──────────────────────────────────────
        md.append("#### Gene model structure\n")
        md.append(img(pd, f"{aid}_gene_length_dist.png"))
        md.append(img(pd, f"{aid}_exons_per_tx.png"))
        md.append(img(pd, f"{aid}_genes_per_scaffold.png"))

        # ── Intron QC ─────────────────────────────────────────────────
        md.append("#### Intron QC\n")
        md.append(img(pd, f"{aid}_intron_qc.png"))

        # ── Repeat annotation ─────────────────────────────────────────
        md.append("#### Repeat annotation\n")
        tbl_file = find_file(args.repeat_tbls, aid, [".tbl"])
        if tbl_file:
            md.append("```\n" + read_text(tbl_file).strip() + "\n```\n")

        tbl_full_file = find_file(args.repeat_tbls_full, aid,
                                   [".buildSummary.tbl"])
        if tbl_full_file:
            md.append("<details>\n<summary>Detailed repeat families (buildSummary.pl)</summary>\n\n")
            md.append("```\n" + read_text(tbl_full_file).strip() + "\n```\n")
            md.append("</details>\n")

        # ── ncRNA ─────────────────────────────────────────────────────
        md.append("#### Non-coding RNA\n")
        md.append(img(pd, f"{aid}_trna_isotypes.png"))
        md.append(img(pd, f"{aid}_ncrna_rfam.png"))

        # ── Telomeres ─────────────────────────────────────────────────
        md.append("#### Telomeres\n")
        md.append(img(pd, f"{aid}_telomere_positions.png"))
        md.append(img(pd, f"{aid}_telomere_density.png"))

        # ── BUSCO ─────────────────────────────────────────────────────
        md.append("#### BUSCO\n")
        busco_file = find_file(args.busco_txts, aid, ["_busco.txt"])
        if busco_file:
            md.append("```\n" + read_text(busco_file).strip() + "\n```\n")

        # ── OMArk ─────────────────────────────────────────────────────
        md.append("#### OMArk\n")
        omark_file = find_file(args.omark_summaries, aid,
                                [".summary.txt"])
        if omark_file:
            om = parse_omark_summary(read_text(omark_file))
            comp = om.get("comp_line", "")
            cons = om.get("cons_line", "")
            clade = om.get("clade", "—")
            n_hogs = om.get("n_hogs", "—")
            n_prot = om.get("n_proteins", "—")

            if comp:
                md.append(f"**Completeness** (clade: {clade}, {n_hogs} conserved HOGs):\n")
                md.append(f"S:{extract_omark_pct(comp, 'S:')}%, "
                          f"D:{extract_omark_pct(comp, 'D:')}% "
                          f"[U:{extract_omark_pct(comp, 'U:')}%, "
                          f"E:{extract_omark_pct(comp, 'E:')}%], "
                          f"M:{extract_omark_pct(comp, 'M:')}%\n\n")
            if cons:
                md.append(f"**Consistency** ({n_prot} proteins):\n")
                md.append(f"A:{extract_omark_pct(cons, 'A:')}% "
                          f"[P:{extract_omark_pct(cons, 'P:')}%, "
                          f"F:{extract_omark_pct(cons, 'F:')}%], "
                          f"I:{extract_omark_pct(cons, 'I:')}%, "
                          f"C:{extract_omark_pct(cons, 'C:')}%, "
                          f"U:{extract_omark_pct(cons, 'U:')}%\n\n")

        omark_det_file = find_file(args.omark_detaileds, aid,
                                    ["_detailed.txt"])
        # Detailed output omitted for brevity — available in qc/omark/ output dir

        # ── AGAT stats ────────────────────────────────────────────────
        agat_file = find_file(args.agat_stats, aid,
                               ["_annotation_stats"])
        if agat_file:
            md.append("<details>\n<summary>AGAT structural statistics</summary>\n\n")
            md.append("```\n" + read_text(agat_file).strip() + "\n```\n")
            md.append("</details>\n")

        # ── Functional annotation ─────────────────────────────────────
        md.append("#### Functional annotation\n")

        n_prot = row.get("n_proteins", "0")
        func_fields = [
            ("pct_ipr", "InterProScan"),
            ("pct_eggnog", "eggNOG"),
            ("pct_diamond", "DIAMOND"),
            ("pct_kofamscan", "KofamScan"),
            ("pct_deepsig", "DeepSig"),
        ]
        md.append(f"| Tool | Coverage |\n|:---|---:|\n")
        for key, label in func_fields:
            md.append(f"| {label} | {row.get(key, '—')}% |\n")
        md.append(f"| **Total proteins** | **{fmt(n_prot)}** |\n\n")

        md.append(img(pd, f"{aid}_functional_coverage.png"))
        md.append(img(pd, f"{aid}_top_products.png"))
        md.append(img(pd, f"{aid}_top_go.png"))
        md.append(img(pd, f"{aid}_go_treemap.png"))
        md.append(img(pd, f"{aid}_go_ontology_composition.png"))
        md.append(img(pd, f"{aid}_top_kegg.png"))
        md.append(img(pd, f"{aid}_top_ipr.png"))

    # ═════════════════════════════════════════════════════════════════════
    # Tuning Recommendations
    # ═════════════════════════════════════════════════════════════════════
    tuning_num = section_num + 1
    md.append(f"---\n## {tuning_num}. Tuning Recommendations\n")

    for row in rows:
        aid = row["assembly_id"]
        flags = []

        # High hypothetical rate
        try:
            pct_hypo = float(row.get("pct_hypothetical", 0))
            if pct_hypo > 50:
                flags.append(f"High hypothetical protein rate ({pct_hypo:.1f}%). "
                             "Consider RNA-seq evidence for BRAKER3 or broader protein databases.")
        except ValueError:
            pass

        # Low BUSCO
        try:
            busco_c = float(row.get("busco_complete", 0))
            if busco_c < 80:
                flags.append(f"Low BUSCO completeness ({busco_c:.1f}%). "
                             "Check repeat masking quality, assembly contiguity, "
                             "or increase BRAKER3 --rounds.")
        except ValueError:
            pass

        # OMArk high unknown
        try:
            omark_u = float(row.get("omark_unknown", 0))
            if omark_u > 25:
                flags.append(f"High OMArk unknown rate ({omark_u:.1f}%). "
                             "Possible over-prediction by BRAKER3. "
                             "Consider filtering with TSEBRA.")
        except ValueError:
            pass

        # OMArk high fragmented
        try:
            omark_f = float(row.get("omark_fragmented", 0))
            if omark_f > 10:
                flags.append(f"High OMArk fragmented rate ({omark_f:.1f}%). "
                             "Many genes shorter than expected — check assembly gaps.")
        except ValueError:
            pass

        # OMArk unexpected duplication
        try:
            omark_du = float(row.get("omark_dup_unexpected", 0))
            if omark_du > 15:
                flags.append(f"High unexpected duplication ({omark_du:.1f}%). "
                             "Possible residual haplotig duplication. "
                             "Cross-reference with BUSCO D% and Hi-C maps.")
        except ValueError:
            pass

        # Short introns
        try:
            si = int(float(row.get("short_introns_lt10", 0)))
            if si > 50:
                flags.append(f"{si:,} introns < 10 bp (likely artifacts). "
                             "Consider post-filtering with AGAT or TSEBRA.")
        except ValueError:
            pass

        # Low repeat masking
        try:
            rmp = float(row.get("repeat_masked_pct", 0))
            if rmp < 20:
                flags.append(f"Low repeat masking ({rmp:.1f}%). "
                             "Fish genomes typically 30-50%. "
                             "Consider more RepeatModeler rounds or slower RM speed.")
        except ValueError:
            pass

        # Low DIAMOND coverage
        try:
            dia_pct = float(row.get("pct_diamond", 0))
            if dia_pct < 30:
                flags.append(f"Low DIAMOND/SwissProt coverage ({dia_pct:.1f}%). "
                             "Consider supplementing with TrEMBL or NR.")
        except ValueError:
            pass

        if flags:
            md.append(f"### `{aid}`\n")
            for f in flags:
                md.append(f"- **{f}**\n")
            md.append("\n")

    if not any(True for row in rows for _ in [1]):  # at least one row exists
        md.append("No major issues detected across assemblies.\n")

    md.append("---\n*Report generated by the genome annotation pipeline.*\n")

    with open(args.outfile, "w") as out:
        out.write("\n".join(md))

    print(f"Report: {n_asm} assemblies → {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()