#!/usr/bin/env python3
"""
generate_annotation_report.py — Unified multi-assembly Markdown annotation report.
"""
import argparse, csv, os, re, sys
from datetime import datetime


# ─── Helpers ─────────────────────────────────────────────────────────────────

def img(pd, md_pd, filename):
    """Markdown image ref. Check existence at pd, write path using md_pd."""
    if os.path.exists(os.path.join(pd, filename)):
        return f"![{filename}]({md_pd}/{filename})"
    return ""

def fmt(n):
    try: return f"{int(float(n)):,}"
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

def extract_omark_pct(text, key):
    m = re.search(re.escape(key) + r"([\d.]+)%", text)
    return m.group(1) if m else "—"

def parse_omark_summary(text):
    r = {}
    for line in text.splitlines():
        ls = line.strip().lstrip("#").strip()
        ll = ls.lower()
        if ls.startswith("S:") and "%" in ls: r["comp_line"] = ls
        if ls.startswith("A:") and "%" in ls: r["cons_line"] = ls
        if "selected clade was" in ll: r["clade"] = ls.split("was")[-1].strip()
        if "conserved hogs" in ll:
            m = re.search(r"(\d[\d,]*)", ls)
            if m: r["n_hogs"] = m.group(1)
        if "there are" in ll and "proteins" in ll:
            m = re.search(r"(\d[\d,]*)\s+proteins", ls)
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
                r["masked_bp"] = m.group(1).replace(",", "")
                r["masked_pct"] = m.group(2)
    return r


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
        vals = [fmt(met[a].get(key, "—")) if is_count else met[a].get(key, "—") for a in aids]
        data.append([label] + vals)
    return table(header, data)

def plot_grid(aids, pd, md_pd, label, tmpl):
    header = [""] + [f"`{a}`" for a in aids]
    imgs = [img(pd, md_pd, tmpl.format(aid=a) + ".png") or "—" for a in aids]
    lines = [
        "| " + " | ".join(header) + " |",
        "|:---|" + "|".join("---" for _ in aids) + "|",
        "| " + label + " | " + " | ".join(imgs) + " |",
    ]
    return "\n".join(lines) + "\n"

def plot_grid_multi(aids, pd, md_pd, rows):
    header = [""] + [f"`{a}`" for a in aids]
    lines = [
        "| " + " | ".join(header) + " |",
        "|:---|" + "|".join("---" for _ in aids) + "|",
    ]
    for label, tmpl in rows:
        imgs = [img(pd, md_pd, tmpl.format(aid=a) + ".png") or "—" for a in aids]
        lines.append("| " + label + " | " + " | ".join(imgs) + " |")
    return "\n".join(lines) + "\n"


# ─── Main ────────────────────────────────────────────────────────────────────

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
    parser.add_argument("--plots_md_path", default=None)
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    pd = args.plots_dir
    md_pd = args.plots_md_path if args.plots_md_path else pd

    met = {}
    with open(args.metrics) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            met[row["assembly_id"]] = row
    aids = sorted(met.keys())

    W = []

    # ═══════════════════════════════════════════════════════════════════════
    # Header
    # ═══════════════════════════════════════════════════════════════════════
    W.append("# Genome Annotation Report\n")
    W.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in aids)}  ")
    W.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    W.append("---\n")

    # ═══════════════════════════════════════════════════════════════════════
    # Annotation outputs
    # ═══════════════════════════════════════════════════════════════════════
    W.append("## Annotation Outputs\n")
    out_files = [
        ("Annotated GFF3",         "{aid}.annotated.gff3"),
        ("Protein sequences",      "{aid}.proteins.fa"),
        ("CDS sequences",          "{aid}.cds.fa"),
        ("Transcript sequences",   "{aid}.transcripts.fa"),
        ("tRNA sequences",         "{aid}.tRNA.fa"),
        ("Functional annotations", "{aid}.functional_annotations.tsv"),
        ("Gene product summary",   "{aid}.gene_product_summary.tsv"),
    ]
    oh = ["File"] + [f"`{a}`" for a in aids]
    od = [[label] + [link(t.format(aid=a)) for a in aids] for label, t in out_files]
    W.append(table(oh, od))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # TOC
    # ═══════════════════════════════════════════════════════════════════════
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
        ("Gene length distribution", "{aid}_gene_length_dist"),
        ("Exons per transcript",     "{aid}_exons_per_tx"),
        ("Intron length QC",         "{aid}_intron_qc"),
        ("Genes per scaffold",       "{aid}_genes_per_scaffold"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 3. Repeat Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 3. Repeat Annotation\n")

    rh = ["Metric"] + [f"`{a}`" for a in aids]
    rd = {}
    for aid in aids:
        tf = find_file(args.repeat_tbls, aid, [".tbl"])
        rd[aid] = parse_repeat_tbl(read_text(tf)) if tf else {}

    rr = []
    for key, label in [("n_sequences","Sequences"), ("total_length","Total length (bp)"),
                        ("gc_pct","GC (%)"), ("masked_bp","Bases masked (bp)"),
                        ("masked_pct","Bases masked (%)")]:
        vals = []
        for a in aids:
            v = rd.get(a, {}).get(key, "—")
            vals.append(fmt(v) if key in ("total_length","masked_bp","n_sequences") else v)
        rr.append([label] + vals)

    # Links row
    lv = []
    for a in aids:
        lv.append(", ".join([
            link(f"repeats/{a}.tbl", "summary"),
            link(f"repeats/{a}.buildSummary.tbl", "families"),
            link(f"repeats/{a}.repeat_landscape.html", "landscape"),
        ]))
    rr.append(["Reports"] + lv)
    W.append(table(rh, rr))
    W.append("")

    W.append(plot_grid(aids, pd, md_pd, "Per-scaffold masking", "{aid}_repeat_masking_per_scaffold"))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 4. Non-coding RNA
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 4. Non-coding RNA\n")
    nh = ["Metric"] + [f"`{a}`" for a in aids]
    nr = [
        ["tRNA genes"] + [fmt(met[a].get("n_trna", "—")) for a in aids],
        ["ncRNA features (Rfam)"] + [fmt(met[a].get("n_ncrna", "—")) for a in aids],
        ["Reports"] + [link("ncRNA/", "ncRNA output") for _ in aids],
    ]
    W.append(table(nh, nr))
    W.append("")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("tRNA isotypes", "{aid}_trna_isotypes"),
        ("Rfam families", "{aid}_ncrna_rfam"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 5. Telomeres & Centromeres
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 5. Telomeres & Centromeres\n")
    W.append(metrics_table(aids, met, [
        ("n_telomere",   "Telomeric regions",  True),
        ("n_centromere", "Centromere features", True),
    ]))
    W.append("")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Telomere positions", "{aid}_telomere_positions"),
        ("Telomere density",   "{aid}_telomere_density"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 6. Quality Assessment
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 6. Quality Assessment\n")

    # BUSCO
    W.append("### BUSCO\n")
    W.append(metrics_table(aids, met, [
        ("busco_complete",   "Complete (%)",    False),
        ("busco_single",     "— Single (%)",   False),
        ("busco_duplicated", "— Duplicated (%)",False),
        ("busco_fragmented", "Fragmented (%)",  False),
        ("busco_missing",    "Missing (%)",     False),
        ("busco_total",      "Total BUSCOs",    True),
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

    # OMArk
    W.append("\n### OMArk\n")
    oi = {}
    for aid in aids:
        of = find_file(args.omark_summaries, aid, [".summary.txt"])
        if of: oi[aid] = parse_omark_summary(read_text(of))

    if oi:
        first = next(iter(oi.values()))
        W.append(f"**Reference clade:** {first.get('clade', '—')}\n")

    W.append("\n#### Completeness (conserved genes)\n")
    h = ["Metric"] + [f"`{a}`" for a in aids]
    cr = [["Conserved HOGs"] + [fmt(oi.get(a, {}).get("n_hogs", "—")) for a in aids]]
    for key, label in [("omark_single","Single (%)"), ("omark_duplicated","Duplicated total (%)"),
                        ("omark_dup_unexpected","— Unexpected (%)"), ("omark_dup_expected","— Expected (%)"),
                        ("omark_missing","Missing (%)")]:
        cr.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(h, cr))
    W.append(f"\n{img(pd, md_pd, 'C_omark_completeness.png')}\n")

    W.append("#### Consistency (proteome quality)\n")
    cr2 = [["Total proteins"] + [fmt(oi.get(a, {}).get("n_proteins", "—")) for a in aids]]
    for key, label in [("omark_consistent","Consistent (%)"), ("omark_partial","— Partial hits (%)"),
                        ("omark_fragmented","— Fragmented (%)"), ("omark_inconsistent","Inconsistent (%)"),
                        ("omark_contaminant","Contaminant (%)"), ("omark_unknown","Unknown (%)")]:
        cr2.append([label] + [met[a].get(key, "—") for a in aids])
    W.append(table(h, cr2))
    W.append(f"\n{img(pd, md_pd, 'C_omark_consistency.png')}\n")
    W.append(f"> Full reports: {link('qc/omark/', 'OMArk output')}\n")

    # AGAT
    W.append("\n### Structural statistics (AGAT)\n")
    kp = [("number of gene","Genes"), ("number of mrna","mRNAs"),
          ("number of exon","Exons"), ("number of cds","CDS features"),
          ("number of intron","Introns"), ("number of single exon gene","Single-exon genes"),
          ("mean gene length","Mean gene length (bp)"), ("mean mrna length","Mean mRNA length (bp)"),
          ("mean cds length","Mean CDS length (bp)"), ("mean exon length","Mean exon length (bp)"),
          ("mean intron length","Mean intron length (bp)"), ("mean exons per","Mean exons/mRNA"),
          ("mean introns per","Mean introns/mRNA"), ("longest gene","Longest gene (bp)"),
          ("shortest gene","Shortest gene (bp)"), ("total gene length","Total gene space (bp)"),
          ("% of genome","% genome in genes")]

    ad = {}
    al = []
    for aid in aids:
        af = find_file(args.agat_stats, aid, ["_annotation_stats"])
        if not af: continue
        ad[aid] = {}
        for line in read_text(af).splitlines():
            ll = line.lower().strip()
            for pat, label in kp:
                if pat in ll:
                    nums = re.findall(r"[\d,.]+", line)
                    if nums:
                        ad[aid][label] = nums[-1]
                        if label not in al: al.append(label)

    if ad:
        ah = ["Metric"] + [f"`{a}`" for a in aids]
        ar = []
        for label in al:
            vals = []
            for a in aids:
                v = ad.get(a, {}).get(label, "—")
                try:
                    if "." not in v and "%" not in label.lower():
                        vals.append(fmt(v))
                    else:
                        vals.append(v)
                except: vals.append(v)
            ar.append([label] + vals)
        W.append(table(ah, ar))

    W.append(f"\n> Full reports: {link('qc/', 'QC output')}\n")

    # ═══════════════════════════════════════════════════════════════════════
    # 7. Functional Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("\n---\n## 7. Functional Annotation\n")
    W.append("### Tool coverage\n")
    W.append(metrics_table(aids, met, [
        ("n_proteins",      "Total proteins",    True),
        ("n_ipr",           "InterProScan hits", True),
        ("pct_ipr",         "InterProScan (%)",  False),
        ("n_eggnog",        "eggNOG hits",       True),
        ("pct_eggnog",      "eggNOG (%)",        False),
        ("n_diamond",       "DIAMOND hits",      True),
        ("pct_diamond",     "DIAMOND (%)",       False),
        ("n_kofamscan",     "KofamScan hits",    True),
        ("pct_kofamscan",   "KofamScan (%)",     False),
        ("n_deepsig",       "DeepSig hits",      True),
        ("pct_deepsig",     "DeepSig (%)",       False),
        ("n_named",         "Named products",    True),
        ("n_hypothetical",  "Hypothetical",      True),
        ("pct_named",       "Named (%)",         False),
        ("pct_hypothetical","Hypothetical (%)",  False),
        ("unique_go_terms", "Unique GO terms",   True),
        ("unique_kegg_kos", "Unique KEGG KOs",   True),
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
            ("pct_hypothetical",    50, None, "High hypothetical rate ({v:.1f}%) — consider RNA-seq evidence or broader protein DBs."),
            ("busco_complete",      None, 80, "Low BUSCO completeness ({v:.1f}%) — check repeat masking, contiguity, or BRAKER3 rounds."),
            ("omark_unknown",       25, None, "High OMArk unknown rate ({v:.1f}%) — possible over-prediction. Consider TSEBRA."),
            ("omark_fragmented",    10, None, "High OMArk fragmented rate ({v:.1f}%) — many genes shorter than expected."),
            ("omark_dup_unexpected", 15, None, "High unexpected duplication ({v:.1f}%) — possible haplotig duplication."),
            ("repeat_masked_pct",   None, 20, "Low repeat masking ({v:.1f}%) — fish genomes typically 30–50%."),
            ("pct_diamond",         None, 30, "Low SwissProt coverage ({v:.1f}%) — consider TrEMBL/NR."),
        ]
        for key, above, below, msg in checks:
            try:
                v = float(row.get(key, 0))
                if (above is not None and v > above) or (below is not None and v < below):
                    flags.append(msg.format(v=v))
            except (ValueError, TypeError): pass
        try:
            si = int(float(row.get("short_introns_lt10", 0)))
            if si > 50: flags.append(f"{si:,} introns < 10 bp — likely gene model artifacts.")
        except (ValueError, TypeError): pass

        if flags:
            any_flags = True
            W.append(f"### `{aid}`\n")
            for f in flags: W.append(f"- {f}")
            W.append("")

    if not any_flags:
        W.append("No major issues detected.\n")

    W.append("---\n*Report generated by the genome annotation pipeline.*\n")

    with open(args.outfile, "w") as out:
        out.write("\n".join(W))
    print(f"Report: {len(aids)} assemblies → {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()