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


def img(pd, md_pd, filename):
    """Return markdown image ref. Check existence at pd, write path as md_pd."""
    path = os.path.join(pd, filename)
    if os.path.exists(path):
        return f"![{filename}]({md_pd}/{filename})"
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
        ll = ls.lower()
        if ls.startswith("S:") and "%" in ls:
            result["comp_line"] = ls
        if ls.startswith("A:") and "%" in ls:
            result["cons_line"] = ls
        if "selected clade was" in ll:
            result["clade"] = ls.split("was")[-1].strip()
        if "conserved hogs" in ll:
            m = re.search(r"(\d[\d,]*)", ls)
            if m:
                result["n_hogs"] = m.group(1)
        if "there are" in ll and "proteins" in ll:
            m = re.search(r"(\d[\d,]*)\s+proteins", ls)
            if m:
                result["n_proteins"] = m.group(1)
    return result


def link(path, text=None):
    """Markdown link to a results path."""
    if text is None:
        text = os.path.basename(path)
    return f"[{text}]({path})"


# ─── Table builders ──────────────────────────────────────────────────────────

def table(header_row, data_rows):
    lines = []
    lines.append("| " + " | ".join(header_row) + " |")
    lines.append("|:---|" + "|".join("---:" for _ in header_row[1:]) + "|")
    for row in data_rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines) + "\n"


def metrics_table(aids, met, fields):
    header = ["Metric"] + [f"`{a}`" for a in aids]
    data = []
    for key, label, is_count in fields:
        vals = []
        for a in aids:
            v = met[a].get(key, "—")
            vals.append(fmt(v) if is_count else v)
        data.append([label] + vals)
    return table(header, data)


def plot_grid(aids, pd, label, filename_template):
    header = [""] + [f"`{a}`" for a in aids]
    imgs = []
    for a in aids:
        fn = filename_template.format(aid=a) + ".png"
        i = img(pd, fn)
        imgs.append(i if i else "—")
    lines = []
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|:---|" + "|".join("---" for _ in aids) + "|")
    lines.append("| " + label + " | " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"


def plot_grid_multi(aids, pd, rows):
    header = [""] + [f"`{a}`" for a in aids]
    lines = []
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|:---|" + "|".join("---" for _ in aids) + "|")
    for label, tmpl in rows:
        imgs = []
        for a in aids:
            fn = tmpl.format(aid=a) + ".png"
            i = img(pd, fn)
            imgs.append(i if i else "—")
        lines.append("| " + label + " | " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"


# ─── Repeat .tbl parser ─────────────────────────────────────────────────────

def parse_repeat_tbl(text):
    """Extract key stats from RM .tbl header."""
    r = {}
    for line in text.splitlines():
        if "sequences:" in line:
            m = re.search(r"(\d+)", line)
            if m: r["n_sequences"] = m.group(1)
        if "total length:" in line:
            m = re.search(r"([\d,]+)\s*bp", line)
            if m: r["total_length"] = m.group(1).replace(",", "")
        if "GC level:" in line:
            m = re.search(r"([\d.]+)\s*%", line)
            if m: r["gc_pct"] = m.group(1)
        if "bases masked:" in line:
            m = re.search(r"([\d,]+)\s*bp\s*\(\s*([\d.]+)\s*%", line)
            if m:
                r["masked_bp"] = m.group(1).replace(",", "")
                r["masked_pct"] = m.group(2)
    return r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--assembly_stats", nargs="+", required=True)
    parser.add_argument("--repeat_tbls", nargs="+", required=True)
    parser.add_argument("--busco_txts", nargs="+", required=True)
    parser.add_argument("--omark_summaries", nargs="+", required=True)
    parser.add_argument("--omark_detaileds", nargs="+", required=True)
    parser.add_argument("--agat_stats", nargs="+", required=True)
    parser.add_argument("--plots_dir", required=True,
                        help="Disk path to plots (for existence checks)")
    parser.add_argument("--plots_md_path", default=None,
                        help="Relative path for markdown image refs (default: same as --plots_dir)")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    pd = args.plots_dir
    md_pd = args.plots_md_path if args.plots_md_path else pd

    met = {}
    with open(args.metrics) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            met[row["assembly_id"]] = row
    aids = sorted(met.keys())
    n_asm = len(aids)

    W = []

    W.append("# Genome Annotation Report\n")
    W.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in aids)}  ")
    W.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    W.append("---\n")

    # ═══════════════════════════════════════════════════════════════════════
    # Annotation outputs table
    # ═══════════════════════════════════════════════════════════════════════
    W.append("## Annotation Outputs\n")
    out_header = ["File"] + [f"`{a}`" for a in aids]
    out_rows = [
        ("Annotated GFF3",        "{aid}.annotated.gff3"),
        ("Protein sequences",     "{aid}.proteins.fa"),
        ("CDS sequences",         "{aid}.cds.fa"),
        ("Transcript sequences",  "{aid}.transcripts.fa"),
        ("tRNA sequences",        "{aid}.tRNA.fa"),
        ("Functional annotations","{aid}.functional_annotations.tsv"),
        ("Gene product summary",  "{aid}.gene_product_summary.tsv"),
    ]
    out_data = []
    for label, path_tmpl in out_rows:
        row = [label]
        for a in aids:
            p = path_tmpl.format(aid=a)
            row.append(link(p))
        out_data.append(row)
    W.append(table(out_header, out_data))
    W.append("")

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
        ("n_genes",              "Protein-coding genes",  True),
        ("n_mrnas",              "mRNA transcripts",      True),
        ("n_proteins",           "Predicted proteins",    True),
        ("n_trna",               "tRNA genes",            True),
        ("n_ncrna",              "ncRNA features",        True),
        ("n_telomere",           "Telomeric regions",     True),
        ("n_centromere",         "Centromere features",   True),
        ("n_repeat",             "Repeat features",       True),
        ("scaffolds_with_genes", "Scaffolds with genes",  True),
        ("repeat_masked_pct",    "Genome masked (%)",     False),
        ("busco_complete",       "BUSCO complete (%)",    False),
        ("omark_consistent",     "OMArk consistent (%)",  False),
        ("n_named",              "Named products",        True),
        ("n_hypothetical",       "Hypothetical proteins", True),
        ("pct_named",            "Named (%)",             False),
        ("unique_go_terms",      "Unique GO terms",       True),
        ("unique_kegg_kos",      "Unique KEGG KOs",       True),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 2. Gene Model Structure
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 2. Gene Model Structure\n")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Gene length distribution",  "{aid}_gene_length_dist"),
        ("Exons per transcript",      "{aid}_exons_per_tx"),
        ("Intron length QC",          "{aid}_intron_qc"),
        ("Genes per scaffold",        "{aid}_genes_per_scaffold"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 3. Repeat Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 3. Repeat Annotation\n")

    # Build repeat summary table with parsed .tbl info
    rep_header = ["Metric"] + [f"`{a}`" for a in aids]
    rep_rows = []
    rep_data = {}
    for aid in aids:
        tbl_file = find_file(args.repeat_tbls, aid, [".tbl"])
        if tbl_file:
            rep_data[aid] = parse_repeat_tbl(read_text(tbl_file))
        else:
            rep_data[aid] = {}

    rep_fields = [
        ("n_sequences",    "Sequences"),
        ("total_length",   "Total length (bp)"),
        ("gc_pct",         "GC (%)"),
        ("masked_bp",      "Bases masked (bp)"),
        ("masked_pct",     "Bases masked (%)"),
    ]
    for key, label in rep_fields:
        vals = []
        for a in aids:
            v = rep_data.get(a, {}).get(key, "—")
            if key in ("total_length", "masked_bp", "n_sequences"):
                vals.append(fmt(v))
            else:
                vals.append(v)
        rep_rows.append([label] + vals)

    # Add links row
    link_vals = []
    for a in aids:
        links = ", ".join([
            link(f"repeats/{a}.tbl", "summary"),
            link(f"repeats/{a}.buildSummary.tbl", "families"),
            link(f"repeats/{a}.repeat_landscape.html", "landscape"),
        ])
        link_vals.append(links)
    rep_rows.append(["Reports"] + link_vals)

    W.append(table(rep_header, rep_rows))
    W.append("")

    W.append(plot_grid(aids, pd, md_pd, "Per-scaffold masking", "{aid}_repeat_masking_per_scaffold"))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 4. Non-coding RNA
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 4. Non-coding RNA\n")
    nc_header = ["Metric"] + [f"`{a}`" for a in aids]
    nc_rows = []
    nc_rows.append(["tRNA genes"] + [fmt(met[a].get("n_trna", "—")) for a in aids])
    nc_rows.append(["ncRNA features (Rfam)"] + [fmt(met[a].get("n_ncrna", "—")) for a in aids])
    nc_rows.append(["Reports"] + [link("ncRNA/", "ncRNA output") for _ in aids])
    W.append(table(nc_header, nc_rows))
    W.append("")

    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("tRNA isotypes",    "{aid}_trna_isotypes"),
        ("Rfam families",    "{aid}_ncrna_rfam"),
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

    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Telomere positions",     "{aid}_telomere_positions"),
        ("Telomere density",       "{aid}_telomere_density"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 6. Quality Assessment
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 6. Quality Assessment\n")

    # ── BUSCO
    W.append("### BUSCO\n")
    W.append(metrics_table(aids, met, [
        ("busco_complete",    "Complete (%)",      False),
        ("busco_single",      "— Single (%)",     False),
        ("busco_duplicated",  "— Duplicated (%)", False),
        ("busco_fragmented",  "Fragmented (%)",    False),
        ("busco_missing",     "Missing (%)",       False),
        ("busco_total",       "Total BUSCOs",      True),
    ]))

    for aid in aids:
        bf = find_file(args.busco_txts, aid, ["_busco.txt"])
        if bf:
            for line in read_text(bf).splitlines():
                if "lineage dataset" in line.lower():
                    W.append(f"\n*{line.strip()}*\n")
                    break
            break

    W.append(f"\n{img(pd, md_pd, 'C_busco.png')}\n")
    W.append(f"> Full reports: {link('qc/busco/', 'BUSCO output')}\n")

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
    header = ["Metric"] + [f"`{a}`" for a in aids]
    rows = []
    rows.append(["Conserved HOGs"] + [
        fmt(omark_info.get(a, {}).get("n_hogs", "—")) for a in aids])
    for key, label in [
        ("omark_single",          "Single (%)"),
        ("omark_duplicated",      "Duplicated total (%)"),
        ("omark_dup_unexpected",  "— Unexpected (%)"),
        ("omark_dup_expected",    "— Expected (%)"),
        ("omark_missing",         "Missing (%)"),
    ]:
        rows.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(header, rows))
    W.append(f"\n{img(pd, md_pd, 'C_omark_completeness.png')}\n")

    W.append("#### Consistency (proteome quality)\n")
    rows2 = []
    rows2.append(["Total proteins"] + [
        fmt(omark_info.get(a, {}).get("n_proteins", "—")) for a in aids])
    for key, label in [
        ("omark_consistent",   "Consistent (%)"),
        ("omark_partial",      "— Partial hits (%)"),
        ("omark_fragmented",   "— Fragmented (%)"),
        ("omark_inconsistent", "Inconsistent (%)"),
        ("omark_contaminant",  "Contaminant (%)"),
        ("omark_unknown",      "Unknown (%)"),
    ]:
        rows2.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(header, rows2))
    W.append(f"\n{img(pd, md_pd, 'C_omark_consistency.png')}\n")
    W.append(f"> Full reports: {link('qc/omark/', 'OMArk output')}\n")

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
            row_vals = []
            for a in aids:
                v = agat_data.get(a, {}).get(label, "—")
                # Add commas to large integers
                try:
                    if "." not in v and "%" not in label.lower():
                        row_vals.append(fmt(v))
                    else:
                        row_vals.append(v)
                except (ValueError, TypeError):
                    row_vals.append(v)
            rows.append([label] + row_vals)
        W.append(table(header, rows))

    W.append(f"\n> Full reports: {link('qc/', 'QC output')}\n")

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
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Tool coverage",           "{aid}_functional_coverage"),
        ("Top products",            "{aid}_top_products"),
        ("Top GO terms",            "{aid}_top_go"),
        ("GO treemap",              "{aid}_go_treemap"),
        ("GO ontology composition", "{aid}_go_ontology_composition"),
        ("Top KEGG orthologs",      "{aid}_top_kegg"),
        ("Top InterPro domains",    "{aid}_top_ipr"),
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
            ("pct_hypothetical", 50, None, "High hypothetical rate ({v:.1f}%) — consider RNA-seq evidence or broader protein DBs."),
            ("busco_complete", None, 80, "Low BUSCO completeness ({v:.1f}%) — check repeat masking, contiguity, or BRAKER3 rounds."),
            ("omark_unknown", 25, None, "High OMArk unknown rate ({v:.1f}%) — possible over-prediction. Consider TSEBRA."),
            ("omark_fragmented", 10, None, "High OMArk fragmented rate ({v:.1f}%) — many genes shorter than expected."),
            ("omark_dup_unexpected", 15, None, "High unexpected duplication ({v:.1f}%) — possible haplotig duplication."),
            ("repeat_masked_pct", None, 20, "Low repeat masking ({v:.1f}%) — fish genomes typically 30–50%."),
            ("pct_diamond", None, 30, "Low SwissProt coverage ({v:.1f}%) — consider TrEMBL/NR."),
        ]

        for key, above, below, msg in checks:
            try:
                v = float(row.get(key, 0))
                trigger = False
                if above is not None and v > above:
                    trigger = True
                if below is not None and v < below:
                    trigger = True
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