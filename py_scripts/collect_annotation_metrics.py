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
    """Count features by source and type from a GFF3 file."""
    counts = defaultdict(int)
    scaffolds_with_genes = set()
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

    # Extract intron lengths for QC
    intron_lengths = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 5 and cols[1] == "BRAKER3" and cols[2] == "intron":
                intron_lengths.append(int(cols[4]) - int(cols[3]) + 1)

    return {
        "n_genes": counts.get(("BRAKER3", "gene"), 0),
        "n_mrnas": counts.get(("BRAKER3", "mRNA"), 0) + counts.get(("BRAKER3", "transcript"), 0),
        "n_exons": counts.get(("BRAKER3", "exon"), 0),
        "n_cds": counts.get(("BRAKER3", "CDS"), 0),
        "n_introns": len(intron_lengths),
        "n_trna": counts.get(("tRNAscan-SE", "gene"), 0) + counts.get(("tRNAscan-SE", "tRNA"), 0),
        "n_ncrna": sum(v for (s, t), v in counts.items() if s == "Infernal"),
        "n_telomere": sum(v for (s, t), v in counts.items() if s == "tidk"),
        "n_centromere": sum(v for (s, t), v in counts.items() if s == "centroAnno"),
        "n_repeat": sum(v for (s, t), v in counts.items() if s == "RepeatMasker"),
        "scaffolds_with_genes": len(scaffolds_with_genes),
        "short_introns_lt10": sum(1 for l in intron_lengths if l < 10),
    }


def parse_gene_summary(path):
    """Count named vs hypothetical products."""
    n_total, n_hypo = 0, 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_total += 1
            if row.get("product", "") == "hypothetical protein":
                n_hypo += 1
    return {"n_genes_total": n_total, "n_hypothetical": n_hypo, "n_named": n_total - n_hypo}


def parse_functional(path):
    """Count functional annotation coverage per tool."""
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
        "n_proteins": n,
        "n_ipr": n_ipr, "n_eggnog": n_egg, "n_diamond": n_dia,
        "n_kofamscan": n_kof, "n_deepsig": n_sig,
        "pct_ipr": round(100 * n_ipr / n, 2) if n else 0,
        "pct_eggnog": round(100 * n_egg / n, 2) if n else 0,
        "pct_diamond": round(100 * n_dia / n, 2) if n else 0,
        "pct_kofamscan": round(100 * n_kof / n, 2) if n else 0,
        "pct_deepsig": round(100 * n_sig / n, 2) if n else 0,
        "unique_go_terms": len(go_terms),
        "unique_kegg_kos": len(kegg_kos),
    }


def parse_repeat_tbl(path):
    """Extract masked % and total masked bp from RM .tbl."""
    text = open(path).read()
    masked_pct = 0.0
    masked_bp = 0
    m = re.search(r"bases masked:\s*([\d,]+)\s*bp\s*\(\s*([\d.]+)\s*%\s*\)", text)
    if m:
        masked_bp = int(m.group(1).replace(",", ""))
        masked_pct = float(m.group(2))
    return {"repeat_masked_bp": masked_bp, "repeat_masked_pct": masked_pct}


def parse_busco(path):
    """Parse BUSCO short summary for C/S/D/F/M percentages."""
    text = open(path).read()
    result = {"busco_complete": 0, "busco_single": 0, "busco_duplicated": 0,
              "busco_fragmented": 0, "busco_missing": 0, "busco_total": 0}
    m = re.search(r"C:([\d.]+)%", text)
    if m:
        result["busco_complete"] = float(m.group(1))
    m = re.search(r"S:([\d.]+)%", text)
    if m:
        result["busco_single"] = float(m.group(1))
    m = re.search(r"D:([\d.]+)%", text)
    if m:
        result["busco_duplicated"] = float(m.group(1))
    m = re.search(r"F:([\d.]+)%", text)
    if m:
        result["busco_fragmented"] = float(m.group(1))
    m = re.search(r"M:([\d.]+)%", text)
    if m:
        result["busco_missing"] = float(m.group(1))
    m = re.search(r"n:\s*(\d+)", text)
    if m:
        result["busco_total"] = int(m.group(1))
    return result


def parse_omark(path):
    """Parse OMArk summary compact format."""
    text = open(path).read()
    result = {"omark_single": 0, "omark_duplicated": 0, "omark_dup_unexpected": 0,
              "omark_dup_expected": 0, "omark_missing": 0,
              "omark_consistent": 0, "omark_inconsistent": 0,
              "omark_contaminant": 0, "omark_unknown": 0,
              "omark_partial": 0, "omark_fragmented": 0}

    def extract(line, key):
        m = re.search(re.escape(key) + r"([\d.]+)%", line)
        return float(m.group(1)) if m else 0.0

    for line in text.splitlines():
        line_s = line.strip().lstrip("#").strip()
        # Completeness pct line
        if line_s.startswith("S:") and "%" in line_s:
            result["omark_single"] = extract(line_s, "S:")
            result["omark_duplicated"] = extract(line_s, "D:")
            result["omark_dup_unexpected"] = extract(line_s, "U:")
            result["omark_dup_expected"] = extract(line_s, "E:")
            result["omark_missing"] = extract(line_s, "M:")
        # Consistency pct line
        if line_s.startswith("A:") and "%" in line_s:
            result["omark_consistent"] = extract(line_s, "A:")
            result["omark_inconsistent"] = extract(line_s, "I:")
            result["omark_contaminant"] = extract(line_s, "C:")
            result["omark_unknown"] = extract(line_s, "U:")
            result["omark_partial"] = extract(line_s, "P:")
            result["omark_fragmented"] = extract(line_s, "F:")

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3s", nargs="+", required=True)
    parser.add_argument("--gene_summaries", nargs="+", required=True)
    parser.add_argument("--functional_tsvs", nargs="+", required=True)
    parser.add_argument("--repeat_tbls", nargs="+", required=True)
    parser.add_argument("--busco_txts", nargs="+", required=True)
    parser.add_argument("--omark_summaries", nargs="+", required=True)
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    # Build per-assembly metrics keyed by assembly_id
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

    # Compute derived metrics
    for aid, m in assemblies.items():
        n_genes = m.get("n_genes_total", m.get("n_genes", 0))
        n_named = m.get("n_named", 0)
        m["pct_named"] = round(100 * n_named / n_genes, 2) if n_genes else 0
        m["pct_hypothetical"] = round(100 * m.get("n_hypothetical", 0) / n_genes, 2) if n_genes else 0

    # Write TSV
    if not assemblies:
        print("WARNING: no assemblies found", file=sys.stderr)
        open(args.outfile, "w").write("assembly_id\n")
        return

    # Determine column order
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