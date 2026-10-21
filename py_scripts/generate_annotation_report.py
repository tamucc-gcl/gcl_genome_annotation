#!/usr/bin/env python3
"""
generate_annotation_report.py — Unified multi-assembly Markdown annotation report.
"""
import argparse, csv, os, re, sys
from datetime import datetime


# ─── Helpers ─────────────────────────────────────────────────────────────────

def img(pd, md_pd, filename):
    if os.path.exists(os.path.join(pd, filename)):
        return f"![{filename}]({md_pd}/{filename})"
    return ""

def fmt(n):
    try: return f"{int(float(n)):,}"
    except (ValueError, TypeError): return str(n)

def fmtf(n, dp=2):
    try: return f"{float(n):.{dp}f}"
    except (ValueError, TypeError): return str(n)

def read_text(path):
    try:
        with open(path) as f: return f.read()
    except (FileNotFoundError, IOError): return ""

def infer_id(path, suffixes):
    base = os.path.basename(path)
    for s in suffixes:
        if base.endswith(s): return base[:-len(s)]
    return os.path.splitext(base)[0]

def find_file(file_list, aid, suffixes):
    for f in file_list:
        if infer_id(f, suffixes) == aid: return f
    return None

def link(path, text=None):
    if text is None: text = os.path.basename(path)
    return f"[{text}]({path})"

def parse_omark_summary(text):
    r = {}
    for line in text.splitlines():
        ls = line.strip().lstrip("#").strip()
        ll = ls.lower()
        if ls.startswith("S:") and "%" in ls: r["comp_line"] = ls
        if ls.startswith("A:") and "%" in ls: r["cons_line"] = ls
        if "selected clade was" in ll: r["clade"] = ls.split("was")[-1].strip()
        if "conserved hogs" in ll:
            m = re.search(r":\s*(\d[\d,]*)", ls)   # matches "is: 14633"
            if m: r["n_hogs"] = m.group(1)
        if "there are" in ll and "proteins" in ll:
            m = re.search(r"there are (\d[\d,]*) proteins", ls)
            if m: r["n_proteins"] = m.group(1)
    return r

def parse_repeat_tbl(text):
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
                r["masked_bp"]  = m.group(1).replace(",", "")
                r["masked_pct"] = m.group(2)
    return r

def parse_structural_qc(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return {}


# ─── Table builders ──────────────────────────────────────────────────────────

def table(header_row, data_rows):
    lines = ["| " + " | ".join(header_row) + " |",
             "|:---|" + "|".join("---:" for _ in header_row[1:]) + "|"]
    for row in data_rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines) + "\n"

def metrics_table(aids, met, fields):
    header = ["Metric"] + [f"`{a}`" for a in aids]
    data = []
    for key, label, is_count in fields:
        vals = [fmt(met[a].get(key, "—")) if is_count else met[a].get(key, "—") for a in aids]
        data.append([label] + vals)
    return table(header, data)

def plot_grid_multi(aids, pd, md_pd, rows):
    header = [""] + [f"`{a}`" for a in aids]
    lines = ["| " + " | ".join(header) + " |",
             "|:---|" + "|".join("---" for _ in aids) + "|"]
    for label, tmpl in rows:
        imgs = [img(pd, md_pd, tmpl.format(aid=a) + ".png") or "—" for a in aids]
        lines.append("| " + label + " | " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"

def plot_grid(aids, pd, md_pd, label, tmpl):
    return plot_grid_multi(aids, pd, md_pd, [(label, tmpl)])


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    from unified_report import main as report_main
    report_main()


if __name__ == "__main__": main()
