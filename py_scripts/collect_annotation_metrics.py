#!/usr/bin/env python3
"""
collect_annotation_metrics.py
=============================
Extract key annotation metrics from per-assembly outputs and write a
single TSV with one row per assembly. This TSV drives the comparison
plots and report.

Inputs are space-separated lists of files (Nextflow collects them).
Assembly IDs are inferred from filenames.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict


def infer_assembly_id(filename, suffixes):
    """Strip known suffixes to recover assembly_id from filename."""
    base = os.path.basename(filename)
    for s in suffixes:
        if base.endswith(s):
            return base[: -len(s)]
    return os.path.splitext(base)[0]


def parse_gff3_counts(path):
    counts = defaultdict(int)
    scaffolds_with_genes = set()
    n_te_flagged = 0
    intron_lengths = []

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            source, ftype = cols[1], cols[2]
            counts[(source, ftype)] += 1
            if source == "BRAKER3" and ftype == "gene":
                scaffolds_with_genes.add(cols[0])
                if "putative_TE_derived" in cols[8]:
                    n_te_flagged += 1
            if source == "BRAKER3" and ftype == "intron":
                intron_lengths.append(int(cols[4]) - int(cols[3]) + 1)

    return {
        "n_genes":               counts.get(("BRAKER3", "gene"), 0),
        "n_mrnas":               counts.get(("BRAKER3", "mRNA"), 0) + counts.get(("BRAKER3", "transcript"), 0),
        "n_exons":               counts.get(("BRAKER3", "exon"), 0),
        "n_cds":                 counts.get(("BRAKER3", "CDS"), 0),
        "n_introns":             len(intron_lengths),
        "n_trna":                counts.get(("tRNAscan-SE", "gene"), 0) + counts.get(("tRNAscan-SE", "tRNA"), 0),
        "n_ncrna":               sum(v for (s, t), v in counts.items() if s == "Infernal"),
        "n_telomere":            sum(v for (s, t), v in counts.items() if s == "tidk"),
        "n_centromere":          sum(v for (s, t), v in counts.items() if s == "centroAnno"),
        "n_repeat":              sum(v for (s, t), v in counts.items() if s == "RepeatMasker"),
        "scaffolds_with_genes":  len(scaffolds_with_genes),
        "short_introns_lt10":    sum(1 for l in intron_lengths if l < 10),
        "n_te_flagged":          n_te_flagged,
    }


def parse_gene_summary(path):
    n_total, n_hypo = 0, 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_total += 1
            if row.get("product", "") == "hypothetical protein":
                n_hypo += 1
    return {"n_genes_total": n_total, "n_hypothetical": n_hypo, "n_named": n_total - n_hypo}


def parse_functional(path):
    n, n_ipr, n_egg, n_dia, n_kof, n_sig = 0, 0, 0, 0, 0, 0
    go_terms, kegg_kos = set(), set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n += 1
            if row.get("ipr_signatures", "-") != "-":
                n_ipr += 1
            if row.get("eggnog_OGs", "-") != "-":
                n_egg += 1
            if row.get("diamond_hit", "-") != "-":
                n_dia += 1
            if row.get("kofam_KO", "-") != "-":
                n_kof += 1
            if row.get("deepsig_type", "-") not in ("-", "Other", ""):
                n_sig += 1
            for field in ("ipr_go", "eggnog_GO"):
                raw = row.get(field, "-")
                if raw and raw != "-":
                    for t in re.split(r"[;,|]+", raw):
                        t = t.strip()
                        if t.startswith("GO:"):
                            go_terms.add(t)
            for field in ("eggnog_KEGG_ko", "kofam_KO"):
                raw = row.get(field, "-")
                if raw and raw != "-":
                    for t in re.split(r"[;,]+", raw):
                        t = t.strip()
                        if t:
                            kegg_kos.add(t)
    return {
        "n_proteins":      n,
        "n_ipr":           n_ipr,   "n_eggnog":      n_egg,
        "n_diamond":       n_dia,   "n_kofamscan":   n_kof,
        "n_deepsig":       n_sig,
        "pct_ipr":         round(100 * n_ipr / n, 2) if n else 0,
        "pct_eggnog":      round(100 * n_egg / n, 2) if n else 0,
        "pct_diamond":     round(100 * n_dia / n, 2) if n else 0,
        "pct_kofamscan":   round(100 * n_kof / n, 2) if n else 0,
        "pct_deepsig":     round(100 * n_sig / n, 2) if n else 0,
        "unique_go_terms": len(go_terms),
        "unique_kegg_kos": len(kegg_kos),
    }


def parse_repeat_tbl(path):
    text = open(path).read()
    masked_pct, masked_bp = 0.0, 0
    m = re.search(r"bases masked:\s*([\d,]+)\s*bp\s*\(\s*([\d.]+)\s*%\s*\)", text)
    if m:
        masked_bp  = int(m.group(1).replace(",", ""))
        masked_pct = float(m.group(2))
    return {"repeat_masked_bp": masked_bp, "repeat_masked_pct": masked_pct}


def parse_busco(path):
    text = open(path).read()
    result = {"busco_complete": 0, "busco_single": 0, "busco_duplicated": 0,
              "busco_fragmented": 0, "busco_missing": 0, "busco_total": 0}
    for key, pat in [("busco_complete","C:"), ("busco_single","S:"),
                     ("busco_duplicated","D:"), ("busco_fragmented","F:"),
                     ("busco_missing","M:")]:
        m = re.search(re.escape(pat) + r"([\d.]+)%", text)
        if m:
            result[key] = float(m.group(1))
    m = re.search(r"n:\s*(\d+)", text)
    if m:
        result["busco_total"] = int(m.group(1))
    return result


def parse_omark(path):
    text = open(path).read()
    result = {k: 0 for k in [
        "omark_single", "omark_duplicated", "omark_dup_unexpected",
        "omark_dup_expected", "omark_missing", "omark_consistent",
        "omark_inconsistent", "omark_contaminant", "omark_unknown",
        "omark_partial", "omark_fragmented",
    ]}

    def extract(line, key):
        m = re.search(re.escape(key) + r"([\d.]+)%", line)
        return float(m.group(1)) if m else 0.0

    for line in text.splitlines():
        ls = line.strip().lstrip("#").strip()
        if ls.startswith("S:") and "%" in ls:
            result["omark_single"]          = extract(ls, "S:")
            result["omark_duplicated"]       = extract(ls, "D:")
            result["omark_dup_unexpected"]   = extract(ls, "U:")
            result["omark_dup_expected"]     = extract(ls, "E:")
            result["omark_missing"]          = extract(ls, "M:")
        if ls.startswith("A:") and "%" in ls:
            result["omark_consistent"]       = extract(ls, "A:")
            result["omark_inconsistent"]     = extract(ls, "I:")
            result["omark_contaminant"]      = extract(ls, "C:")
            result["omark_unknown"]          = extract(ls, "U:")
            result["omark_partial"]          = extract(ls, "P:")
            result["omark_fragmented"]       = extract(ls, "F:")
    return result


def parse_structural_qc(path):
    """
    Read the single-row structural_qc.tsv produced by STRUCTURAL_QC.
    All column names are already well-named; we pick the subset most
    useful for cross-assembly comparison and report flagging.
    """
    WANTED = {
        "pct_genes_pass", "pct_genes_warn", "pct_genes_fail",
        "n_high_confidence", "pct_high_confidence",
        "pct_canonical_splice_sites", "n_total_noncanonical",
        "pct_single_exon_mrna", "n_short_introns_lt60",
        "pct_short_introns_lt60",
        "pct_mrnas_with_start_m", "pct_mrnas_with_internal_stop",
        "pct_mrnas_no_stop_codon",
        "mean_protein_length_aa", "median_protein_length_aa",
        "mean_exons_per_mrna", "median_exons_per_mrna",
        "omark_class_distribution",
        "mean_aed", "median_aed", "pct_aed_lt03", "pct_aed_lt05", "pct_aed_no_evidence"
    }
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return {k: v for k, v in row.items() if k in WANTED}
    return {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3s",              nargs="+", required=True)
    parser.add_argument("--gene_summaries",     nargs="+", required=True)
    parser.add_argument("--functional_tsvs",    nargs="+", required=True)
    parser.add_argument("--repeat_tbls",        nargs="+", required=True)
    parser.add_argument("--busco_txts",         nargs="+", required=True)
    parser.add_argument("--omark_summaries",    nargs="+", required=True)
    parser.add_argument("--structural_qc_tsvs", nargs="+", required=True)
    parser.add_argument("--outfile",            required=True)
    args = parser.parse_args()

    assemblies = {}

    for path in args.gff3s:
        aid = infer_assembly_id(path, [".annotated.gff3"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_gff3_counts(path))

    for path in args.gene_summaries:
        aid = infer_assembly_id(path, [".gene_product_summary.tsv"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_gene_summary(path))

    for path in args.functional_tsvs:
        aid = infer_assembly_id(path, [".functional_annotations.tsv"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_functional(path))

    for path in args.repeat_tbls:
        aid = infer_assembly_id(path, [".tbl"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_repeat_tbl(path))

    for path in args.busco_txts:
        aid = infer_assembly_id(path, ["_busco.txt"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_busco(path))

    for path in args.omark_summaries:
        aid = infer_assembly_id(path, [".summary.txt"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_omark(path))

    for path in args.structural_qc_tsvs:
        aid = infer_assembly_id(path, [".structural_qc.tsv"])
        assemblies.setdefault(aid, {})["assembly_id"] = aid
        assemblies[aid].update(parse_structural_qc(path))

    # Derived metrics
    for aid, m in assemblies.items():
        n_genes = m.get("n_genes_total", m.get("n_genes", 0))
        n_named = m.get("n_named", 0)
        m["pct_named"]        = round(100 * n_named / n_genes, 2) if n_genes else 0
        m["pct_hypothetical"] = round(100 * m.get("n_hypothetical", 0) / n_genes, 2) if n_genes else 0
        n_te = m.get("n_te_flagged", 0)
        m["pct_te_flagged"]   = round(100 * n_te / n_genes, 2) if n_genes else 0

    if not assemblies:
        print("WARNING: no assemblies found", file=sys.stderr)
        open(args.outfile, "w").write("assembly_id\n")
        return

    all_keys = ["assembly_id"]
    for m in assemblies.values():
        for k in m:
            if k not in all_keys:
                all_keys.append(k)

    with open(args.outfile, "w") as out:
        out.write("\t".join(all_keys) + "\n")
        for aid in sorted(assemblies):
            m = assemblies[aid]
            vals = [str(m.get(k, "")) for k in all_keys]
            out.write("\t".join(vals) + "\n")

    print(f"Collected metrics for {len(assemblies)} assemblies → {args.outfile}",
          file=sys.stderr)


if __name__ == "__main__":
    main()