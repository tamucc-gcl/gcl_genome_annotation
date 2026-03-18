#!/usr/bin/env python3
"""
generate_annotation_report.py
=============================
Unified multi-assembly Markdown annotation report.
"""

import argparse
import csv
import os
import re
import sys
from datetime import datetime


def img(pd, filename):
    path = os.path.join(pd, filename)
    if os.path.exists(path):
        return f"![{filename}]({pd}/{filename})"
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


# ─── Table builders ──────────────────────────────────────────────────────────
# These return a single string with no trailing blank lines between rows.

def table(header_row, data_rows):
    """Build a markdown table from header list and list-of-lists data."""
    lines = []
    lines.append("| " + " | ".join(header_row) + " |")
    lines.append("|:---|" + "|".join("---:" for _ in header_row[1:]) + "|")
    for row in data_rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines) + "\n"


def metrics_table(aids, met, fields):
    """Cross-assembly metrics table: rows=metrics, cols=assemblies."""
    header = ["Metric"] + [f"`{a}`" for a in aids]
    data = []
    for key, label, is_count in fields:
        vals = []
        for a in aids:
            v = met[a].get(key, "—")
            vals.append(fmt(v) if is_count else v)
        data.append([label] + vals)
    return table(header, data)


def plot_grid(aids, pd, filename_template):
    """Render plots in a markdown table with one column per assembly."""
    header = [f"`{a}`" for a in aids]
    imgs = []
    for a in aids:
        i = img(pd, filename_template.format(aid=a))
        imgs.append(i if i else "—")
    lines = []
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "|".join("---" for _ in aids) + "|")
    lines.append("| " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"


def plot_grid_multi(aids, pd, filename_templates):
    """Multiple rows of plots in a table, one column per assembly."""
    header = [f"`{a}`" for a in aids]
    lines = []
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "|".join("---" for _ in aids) + "|")
    for tmpl in filename_templates:
        imgs = []
        for a in aids:
            i = img(pd, tmpl.format(aid=a))
            imgs.append(i if i else "—")
        lines.append("| " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"


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

    W = []  # output lines

    W.append("# Genome Annotation Report\n")
    W.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in aids)}  ")
    W.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    W.append("---\n")

    # ── TOC
    W.append("""## Table of Contents

1. [Overview](#1-overview)
2. [Gene Model Structure](#2-gene-model-structure)
3. [Repeat Annotation](#3-repeat-annotation)
4. [Non-coding RNA](#4-non-coding-rna)
5. [Telomeres & Centromeres](#5-telomeres--centromeres)
6. [Quality Assessment](#6-quality-assessment)
7. [Functional Annotation](#7-functional-annotation)
8. [Tuning Recommendations](#8-tuning-recommendations)
""")

    # ═══════════════════════════════════════════════════════════════════════
    # 1. Overview
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 1. Overview\n")
    W.append(metrics_table(aids, met, [
        ("n_genes",           "Protein-coding genes",  True),
        ("n_mrnas",           "mRNA transcripts",      True),
        ("n_proteins",        "Predicted proteins",    True),
        ("n_trna",            "tRNA genes",            True),
        ("n_ncrna",           "ncRNA features",        True),
        ("n_telomere",        "Telomeric regions",     True),
        ("n_centromere",      "Centromere features",   True),
        ("n_repeat",          "Repeat features",       True),
        ("scaffolds_with_genes", "Scaffolds with genes", True),
        ("repeat_masked_pct", "Genome masked (%)",     False),
        ("busco_complete",    "BUSCO complete (%)",    False),
        ("omark_consistent",  "OMArk consistent (%)",  False),
        ("n_named",           "Named products",        True),
        ("n_hypothetical",    "Hypothetical proteins", True),
        ("pct_named",         "Named (%)",             False),
        ("unique_go_terms",   "Unique GO terms",       True),
        ("unique_kegg_kos",   "Unique KEGG KOs",       True),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 2. Gene Model Structure
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 2. Gene Model Structure\n")
    W.append(plot_grid_multi(aids, pd, [
        "{aid}_gene_length_dist",
        "{aid}_exons_per_tx",
        "{aid}_intron_qc",
        "{aid}_genes_per_scaffold",
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 3. Repeat Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 3. Repeat Annotation\n")
    W.append(metrics_table(aids, met, [
        ("repeat_masked_pct", "Genome masked (%)", False),
        ("repeat_masked_bp",  "Masked bp",         True),
        ("n_repeat",          "Repeat features",   True),
    ]))
    W.append("")

    W.append(plot_grid(aids, pd, "{aid}_repeat_masking_per_scaffold"))
    W.append("")

    # Per-assembly compact summary
    for aid in aids:
        tbl_file = find_file(args.repeat_tbls, aid, [".tbl"])
        if not tbl_file:
            continue
        tbl_text = read_text(tbl_file).strip()
        summary_lines = []
        for line in tbl_text.splitlines():
            summary_lines.append(line)
            if "bases masked" in line.lower():
                break
        if summary_lines:
            W.append(f"**{aid}:**\n```\n" + "\n".join(summary_lines) + "\n```\n")

    W.append(f"> Full repeat tables: `results/repeats/`\n")

    # ═══════════════════════════════════════════════════════════════════════
    # 4. Non-coding RNA
    # ═══════════════════════════════════════════════════════════════════════
    W.append("\n---\n## 4. Non-coding RNA\n")
    W.append(metrics_table(aids, met, [
        ("n_trna",  "tRNA genes",           True),
        ("n_ncrna", "ncRNA features (Rfam)", True),
    ]))
    W.append("")

    W.append(plot_grid_multi(aids, pd, [
        "{aid}_trna_isotypes",
        "{aid}_ncrna_rfam",
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 5. Telomeres & Centromeres
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 5. Telomeres & Centromeres\n")
    W.append(metrics_table(aids, met, [
        ("n_telomere",    "Telomeric regions",   True),
        ("n_centromere",  "Centromere features",  True),
    ]))
    W.append("")

    W.append(plot_grid_multi(aids, pd, [
        "{aid}_telomere_positions",
        "{aid}_telomere_density",
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 6. Quality Assessment
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 6. Quality Assessment\n")

    # ── BUSCO
    W.append("### BUSCO\n")
    W.append(metrics_table(aids, met, [
        ("busco_complete",    "Complete (%)",    False),
        ("busco_single",      "— Single (%)",   False),
        ("busco_duplicated",  "— Duplicated (%)", False),
        ("busco_fragmented",  "Fragmented (%)",  False),
        ("busco_missing",     "Missing (%)",     False),
        ("busco_total",       "Total BUSCOs",    True),
    ]))

    # Lineage
    for aid in aids:
        bf = find_file(args.busco_txts, aid, ["_busco.txt"])
        if bf:
            for line in read_text(bf).splitlines():
                if "lineage dataset" in line.lower():
                    W.append(f"\n*{line.strip()}*\n")
                    break
            break

    W.append(f"\n{img(pd, 'C_busco.png')}\n")
    W.append(f"> Full BUSCO reports: `results/qc/busco/`\n")

    # ── OMArk
    W.append("\n### OMArk\n")

    omark_info = {}
    for aid in aids:
        of = find_file(args.omark_summaries, aid, [".summary.txt"])
        if of:
            omark_info[aid] = parse_omark_summary(read_text(of))

    if omark_info:
        first = next(iter(omark_info.values()))
        W.append(f"**Reference clade:** {first.get('clade', '—')}\n")

    W.append("\n#### Completeness (conserved genes)\n")

    # Build completeness table manually to include n_hogs
    header = ["Metric"] + [f"`{a}`" for a in aids]
    rows = []
    rows.append(["Conserved HOGs"] + [omark_info.get(a, {}).get("n_hogs", "—") for a in aids])
    for key, label in [
        ("omark_single", "Single (%)"),
        ("omark_duplicated", "Duplicated total (%)"),
        ("omark_dup_unexpected", "— Unexpected (%)"),
        ("omark_dup_expected", "— Expected (%)"),
        ("omark_missing", "Missing (%)"),
    ]:
        rows.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(header, rows))
    W.append(f"\n{img(pd, 'C_omark_completeness.png')}\n")

    W.append("#### Consistency (proteome quality)\n")
    rows2 = []
    rows2.append(["Total proteins"] + [omark_info.get(a, {}).get("n_proteins", "—") for a in aids])
    for key, label in [
        ("omark_consistent", "Consistent (%)"),
        ("omark_partial", "— Partial hits (%)"),
        ("omark_fragmented", "— Fragmented (%)"),
        ("omark_inconsistent", "Inconsistent (%)"),
        ("omark_contaminant", "Contaminant (%)"),
        ("omark_unknown", "Unknown (%)"),
    ]:
        rows2.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(header, rows2))
    W.append(f"\n{img(pd, 'C_omark_consistency.png')}\n")
    W.append(f"> Full OMArk reports: `results/qc/omark/`\n")

    # ── AGAT
    W.append("\n### Structural statistics (AGAT)\n")

    key_patterns = [
        ("number of gene",              "Genes"),
        ("number of mrna",              "mRNAs"),
        ("number of exon",              "Exons"),
        ("number of cds",               "CDS features"),
        ("number of intron",            "Introns"),
        ("number of single exon gene",  "Single-exon genes"),
        ("mean gene length",            "Mean gene length (bp)"),
        ("mean mrna length",            "Mean mRNA length (bp)"),
        ("mean cds length",             "Mean CDS length (bp)"),
        ("mean exon length",            "Mean exon length (bp)"),
        ("mean intron length",          "Mean intron length (bp)"),
        ("mean exons per",              "Mean exons/mRNA"),
        ("mean introns per",            "Mean introns/mRNA"),
        ("longest gene",                "Longest gene (bp)"),
        ("shortest gene",               "Shortest gene (bp)"),
        ("total gene length",           "Total gene space (bp)"),
        ("% of genome",                 "% genome in genes"),
    ]

    agat_data = {}
    agat_labels_ordered = []
    for aid in aids:
        agat_file = find_file(args.agat_stats, aid, ["_annotation_stats"])
        if not agat_file:
            continue
        agat_text = read_text(agat_file)
        agat_data[aid] = {}
        for line in agat_text.splitlines():
            ll = line.lower().strip()
            for pattern, label in key_patterns:
                if pattern in ll:
                    nums = re.findall(r"[\d,.]+", line)
                    if nums:
                        agat_data[aid][label] = nums[-1]
                        if label not in agat_labels_ordered:
                            agat_labels_ordered.append(label)

    if agat_data:
        header = ["Metric"] + [f"`{a}`" for a in aids]
        rows = []
        for label in agat_labels_ordered:
            rows.append([label] + [agat_data.get(a, {}).get(label, "—") for a in aids])
        W.append(table(header, rows))

    W.append(f"\n> Full AGAT reports: `results/qc/`\n")

    # ═══════════════════════════════════════════════════════════════════════
    # 7. Functional Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("\n---\n## 7. Functional Annotation\n")

    W.append("### Tool coverage\n")
    W.append(metrics_table(aids, met, [
        ("n_proteins",      "Total proteins",      True),
        ("n_ipr",           "InterProScan hits",   True),
        ("pct_ipr",         "InterProScan (%)",    False),
        ("n_eggnog",        "eggNOG hits",         True),
        ("pct_eggnog",      "eggNOG (%)",          False),
        ("n_diamond",       "DIAMOND hits",        True),
        ("pct_diamond",     "DIAMOND (%)",         False),
        ("n_kofamscan",     "KofamScan hits",      True),
        ("pct_kofamscan",   "KofamScan (%)",       False),
        ("n_deepsig",       "DeepSig hits",        True),
        ("pct_deepsig",     "DeepSig (%)",         False),
        ("n_named",         "Named products",      True),
        ("n_hypothetical",  "Hypothetical",        True),
        ("pct_named",       "Named (%)",           False),
        ("pct_hypothetical","Hypothetical (%)",    False),
        ("unique_go_terms", "Unique GO terms",     True),
        ("unique_kegg_kos", "Unique KEGG KOs",     True),
    ]))
    W.append("")

    W.append("### Details\n")
    W.append(plot_grid_multi(aids, pd, [
        "{aid}_functional_coverage",
        "{aid}_top_products",
        "{aid}_top_go",
        "{aid}_go_treemap",
        "{aid}_go_ontology_composition",
        "{aid}_top_kegg",
        "{aid}_top_ipr",
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 8. Tuning Recommendations
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 8. Tuning Recommendations\n")

    any_flags = False
    for aid in aids:
        row = met[aid]
        flags = []

        checks = [
            ("pct_hypothetical", 50, "High hypothetical rate ({v:.1f}%) — consider RNA-seq evidence or broader protein DBs."),
            ("busco_complete", None, "Low BUSCO completeness ({v:.1f}%) — check repeat masking, contiguity, or BRAKER3 rounds."),
            ("omark_unknown", 25, "High OMArk unknown rate ({v:.1f}%) — possible over-prediction. Consider TSEBRA."),
            ("omark_fragmented", 10, "High OMArk fragmented rate ({v:.1f}%) — many genes shorter than expected."),
            ("omark_dup_unexpected", 15, "High unexpected duplication ({v:.1f}%) — possible haplotig duplication."),
            ("repeat_masked_pct", None, "Low repeat masking ({v:.1f}%) — fish genomes typically 30–50%."),
            ("pct_diamond", None, "Low SwissProt coverage ({v:.1f}%) — consider TrEMBL/NR."),
        ]

        for key, thresh, msg in checks:
            try:
                v = float(row.get(key, 0))
                trigger = False
                if key == "busco_complete":
                    trigger = v < 80
                elif key == "repeat_masked_pct":
                    trigger = v < 20
                elif key == "pct_diamond":
                    trigger = v < 30
                elif thresh is not None:
                    trigger = v > thresh
                if trigger:
                    flags.append(msg.format(v=v))
            except (ValueError, TypeError):
                pass

        try:
            si = int(float(row.get("short_introns_lt10", 0)))
            if si > 50:
                flags.append(f"{si:,} introns < 10 bp — likely gene model artifacts.")
        except (ValueError, TypeError):
            pass

        if flags:
            any_flags = True
            W.append(f"### `{aid}`\n")
            for f in flags:
                W.append(f"- {f}")
            W.append("")

    if not any_flags:
        W.append("No major issues detected.\n")

    W.append("---\n*Report generated by the genome annotation pipeline.*\n")

    with open(args.outfile, "w") as out:
        out.write("\n".join(W))

    print(f"Report: {n_asm} assemblies → {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()