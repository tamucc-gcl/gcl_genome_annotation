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
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics",            required=True)
    parser.add_argument("--assembly_stats",     nargs="+", required=True)
    parser.add_argument("--repeat_tbls",        nargs="+", required=True)
    parser.add_argument("--busco_txts",         nargs="+", required=True)
    parser.add_argument("--omark_summaries",    nargs="+", required=True)
    parser.add_argument("--omark_detaileds",    nargs="+", required=True)
    parser.add_argument("--agat_stats",         nargs="+", required=True)
    parser.add_argument("--structural_qc_tsvs", nargs="+", required=True)
    parser.add_argument("--pipeline_params",    required=True)
    parser.add_argument("--plots_dir",          required=True)
    parser.add_argument("--plots_md_path",      default=None)
    parser.add_argument("--outfile",            required=True)
    args = parser.parse_args()

    pd    = args.plots_dir
    md_pd = args.plots_md_path if args.plots_md_path else pd

    # Load metrics
    met = {}
    with open(args.metrics) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            met[row["assembly_id"]] = row
    aids = sorted(met.keys())

    # Load pipeline params
    params = {}
    with open(args.pipeline_params) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if "\t" in line:
                k, v = line.split("\t", 1)
                params[k.strip()] = v.strip()

    # Load per-assembly structural QC
    sqc = {}
    for path in args.structural_qc_tsvs:
        aid = infer_id(path, [".structural_qc.tsv"])
        sqc[aid] = parse_structural_qc(path)

    W = []

    # ═══════════════════════════════════════════════════════════════════════
    # Header
    # ═══════════════════════════════════════════════════════════════════════
    W.append("# Genome Annotation Report\n")
    W.append(f"**Assemblies:** {', '.join(f'`{a}`' for a in aids)}  ")
    W.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}  ")
    W.append(f"**Species:** {params.get('taxa_name', 'unset')}  ")
    if params.get('ncbi_speciesID', 'none') != 'none':
        W.append(f"**NCBI TaxID:** {params['ncbi_speciesID']}  ")
    W.append("\n---\n")

    # ═══════════════════════════════════════════════════════════════════════
    # Pipeline configuration & databases
    # ═══════════════════════════════════════════════════════════════════════
    W.append("## Pipeline Configuration\n")

    W.append("### Tools & Databases\n")
    db_rows = [
        ("Repeat annotation", [
            ("RepeatModeler2 + RepeatMasker", "de novo repeat families + masking"),
            (f"Dfam v{params.get('dfam_version','?')}", f"Taxon: `{params.get('taxa_name','?')}`"),
            (f"Terrier (threshold {params.get('terrier_threshold','?')})", "Unknown repeat reclassification"),
            (f"DIAMOND filter (pident ≥ {params.get('repeat_filter_pident','?')}%, e-value ≤ 1e-10)",
             "Gene-derived repeat removal vs SwissProt + RepeatPeps"),
        ]),
        ("Structural annotation", [
            (f"BRAKER3 (rounds={params.get('braker_rounds','?')})", "Protein-guided gene prediction"),
            (f"OrthoDB v12 `{params.get('orthodb_clade','none')}`", "Training protein database"),
            (f"Reference proteome: `{params.get('ref_proteome_accession','none')}`",
             "Close-relative protein evidence"),
            ("tRNAscan-SE", "tRNA annotation"),
            (f"Infernal / Rfam `{params.get('rfam_version','?')}`", "ncRNA annotation"),
            ("tidk", f"Telomere detection (motif: `{params.get('telomere_motif','TTAGGG')}`)"),
            ("centroAnno", "Centromere annotation"),
        ]),
        ("Quality assessment", [
            (f"BUSCO lineage: `{params.get('busco_lineage','?')}`", "Gene space completeness"),
            (f"OMArk database: `{params.get('oma_db_name','?')}`",
             f"Proteome consistency (taxid {params.get('ncbi_speciesID','none')})"),
            ("AGAT", "Structural annotation statistics"),
            ("Splice-site audit", "GT-AG / GC-AG / AT-AC dinucleotide check per intron"),
        ]),
        ("Functional annotation", [
            (f"InterProScan `{params.get('interproscan_version','?')}`", "Domain / family annotation"),
            (f"eggNOG-mapper `{params.get('eggnog_version','?')}`", "Ortholog groups + GO + KEGG"),
            (f"DIAMOND blastp vs SwissProt (`{params.get('diamond_version','?')}`)", "Best SwissProt hit"),
            (f"KofamScan (`{params.get('kofamscan_version','?')}`)", "KEGG Ortholog assignment"),
            ("DeepSig", "Signal peptide prediction"),
        ]),
    ]

    for category, entries in db_rows:
        W.append(f"#### {category}\n")
        rows = [[tool, note] for tool, note in entries]
        W.append(table(["Tool / Database", "Notes"], rows))
        W.append("")

    W.append("### Key Filtering Parameters\n")
    param_rows = [
        ("Minimum contig length",       params.get("min_contig_length", "?") + " bp"),
        ("Min CDS length (single-exon)",params.get("filter_min_cds_length", "?") + " bp"),
        ("Min gene span",               params.get("filter_min_gene_length", "?") + " bp"),
        ("Min intron length (filter)",  params.get("filter_min_intron_length", "?") + " bp"),
        ("TE gene overlap threshold",   params.get("te_flag_overlap_frac", "?")),
        ("RepeatMasker chunks",         params.get("n_rm_chunks", "?")),
    ]
    W.append(table(["Parameter", "Value"], param_rows))
    W.append("")

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
        ("High-confidence genes",  "qc/structural/{aid}.high_confidence_genes.txt"),
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
3. [Structural Quality Assessment](#3-structural-quality-assessment)
4. [Repeat Annotation](#4-repeat-annotation)
5. [Non-coding RNA](#5-non-coding-rna)
6. [Telomeres & Centromeres](#6-telomeres--centromeres)
7. [Completeness & Quality](#7-completeness--quality)
8. [Functional Annotation](#8-functional-annotation)
9. [Tuning Recommendations](#9-tuning-recommendations)
""")

    # ═══════════════════════════════════════════════════════════════════════
    # 1. Overview
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 1. Overview\n")
    W.append(metrics_table(aids, met, [
        ("n_genes",              "Protein-coding genes",    True),
        ("n_mrnas",              "mRNA transcripts",        True),
        ("n_proteins",           "Predicted proteins",      True),
        ("n_trna",               "tRNA genes",              True),
        ("n_ncrna",              "ncRNA features",          True),
        ("n_telomere",           "Telomeric regions",       True),
        ("n_centromere",         "Centromere features",     True),
        ("n_repeat",             "Repeat features",         True),
        ("scaffolds_with_genes", "Scaffolds with genes",    True),
        ("repeat_masked_pct",    "Genome masked (%)",       False),
        ("busco_complete",       "BUSCO complete (%)",      False),
        ("omark_consistent",     "OMArk consistent (%)",    False),
        ("n_high_confidence",    "High-confidence genes",   True),
        ("pct_high_confidence",  "High-confidence (%)",     False),
        ("n_named",              "Named products",          True),
        ("n_hypothetical",       "Hypothetical proteins",   True),
        ("pct_named",            "Named (%)",               False),
        ("n_te_flagged",         "TE-flagged genes",        True),
        ("pct_te_flagged",       "TE-flagged (%)",          False),
        ("unique_go_terms",      "Unique GO terms",         True),
        ("unique_kegg_kos",      "Unique KEGG KOs",         True),
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
    # 3. Structural Quality Assessment
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 3. Structural Quality Assessment\n")

    W.append("""\
Per-gene quality flags are assigned by inspecting each predicted transcript for:
internal stop codons, missing start Met, short proteins (<100 aa), short introns (<60 bp),
non-canonical splice sites, and single-exon models without functional support.
The **worst** flag across all transcripts determines the gene-level flag.
**High-confidence** genes must have all transcripts flagged PASS and must not be TE-flagged.
""")

    # Cross-assembly quality flag table
    W.append("### Quality flag summary\n")
    W.append(metrics_table(aids, met, [
        ("pct_genes_pass",  "PASS (%)",  False),
        ("pct_genes_warn",  "WARN (%)",  False),
        ("pct_genes_fail",  "FAIL (%)",  False),
        ("n_high_confidence",   "High-confidence genes",   True),
        ("pct_high_confidence", "High-confidence (%)",     False),
    ]))
    W.append("")

    # Protein sequence integrity
    W.append("### Protein sequence integrity\n")
    W.append(metrics_table(aids, met, [
        ("pct_mrnas_with_start_m",       "Has start Met (%)",       False),
        ("pct_mrnas_no_stop_codon",      "Missing stop codon (%)",  False),
        ("pct_mrnas_with_internal_stop", "Internal stop codon (%)", False),
        ("mean_protein_length_aa",       "Mean protein length (aa)",False),
        ("median_protein_length_aa",     "Median protein length (aa)", False),
    ]))
    W.append("""\
> **Guidance:** Start Met < 95% indicates fragmented assembly or annotation at scaffold edges.
> Internal stop codons > 1% suggest gene model artifacts or assembly errors.
> Missing stop codons are common for genes near scaffold ends — usually benign.
""")

    # Splice site audit
    W.append("### Splice-site audit\n")
    W.append(metrics_table(aids, met, [
        ("pct_canonical_splice_sites", "Canonical splice sites (%)", False),
        ("n_total_noncanonical",       "Non-canonical introns",      True),
        ("n_short_introns_lt60",       "Short introns < 60 bp",      True),
        ("pct_short_introns_lt60",     "Short introns (%)",          False),
    ]))
    W.append("""\
> **Guidance:** Canonical (GT-AG/GC-AG/AT-AC) rate should be >98% in a well-assembled genome.
> Short introns < 60 bp are likely gene model artifacts from GeneMark training.
> Both are strong indicators for filtering the gene set before downstream analysis.
""")

    # Evidence support (AED)
    W.append("### Evidence support (AED)\n")
    W.append("""\
Annotation Edit Distance (AED) measures congruency between each BRAKER gene model
and independently aligned protein evidence (miniprot).  AED = 0 indicates perfect
agreement with evidence; AED = 1 indicates no overlapping evidence.  Values < 0.5
are considered well-supported; < 0.3 high-quality (Eilbeck et al. 2009).
""")
    W.append(metrics_table(aids, met, [
        ("mean_aed",            "Mean AED",                 False),
        ("median_aed",          "Median AED",               False),
        ("pct_aed_lt03",        "AED < 0.3 (%)",            False),
        ("pct_aed_lt05",        "AED < 0.5 (%)",            False),
        ("pct_aed_no_evidence", "No evidence, AED=1 (%)",   False),
    ]))
    W.append("")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Cumulative AED", "{aid}_aed_cumulative"),
        ("AED histogram",  "{aid}_aed_histogram"),
    ]))
    W.append("")

    # Gene structure statistics
    W.append("### Gene structure statistics\n")
    W.append(metrics_table(aids, met, [
        ("mean_exons_per_mrna",    "Mean exons/mRNA",    False),
        ("median_exons_per_mrna",  "Median exons/mRNA",  False),
        ("pct_single_exon_mrna",   "Single-exon mRNAs (%)", False),
        ("n_te_flagged",           "TE-flagged genes",   True),
        ("pct_te_flagged",         "TE-flagged (%)",     False),
    ]))
    W.append("""\
> **Guidance:** Teleost genomes typically have 6–10 exons/transcript and 15–25% single-exon genes.
> Higher single-exon rates often indicate TE contamination of the gene set.
""")

    # Per-assembly detail tables
    W.append("### Per-assembly detail\n")
    for aid in aids:
        q = sqc.get(aid, {})
        if not q:
            continue
        W.append(f"#### `{aid}`\n")

        # OMArk class distribution breakdown
        omark_dist = q.get("omark_class_distribution", "-")
        if omark_dist and omark_dist != "-":
            W.append(f"**OMArk per-protein class distribution:** `{omark_dist}`\n")

        # Splice-site class breakdown from functional TSV is not in sqc directly;
        # the donor_acceptor_summary column in the TSV is per-mRNA. The aggregate
        # is captured in pct_canonical_splice_sites already displayed above.

        W.append(f"> Full gene-level QC table: {link(f'qc/structural/{aid}.structural_qc_genes.tsv', 'per-gene TSV')}\n")
        W.append(f"> High-confidence gene list: {link(f'qc/structural/{aid}.high_confidence_genes.txt', 'gene list')}\n")

    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 4. Repeat Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 4. Repeat Annotation\n")

    rh = ["Metric"] + [f"`{a}`" for a in aids]
    rd = {}
    for aid in aids:
        tf = find_file(args.repeat_tbls, aid, [".tbl"])
        rd[aid] = parse_repeat_tbl(read_text(tf)) if tf else {}

    rr = []
    for key, label in [("n_sequences","Sequences"), ("total_length","Total length (bp)"),
                        ("gc_pct","GC (%)"), ("masked_bp","Bases masked (bp)"),
                        ("masked_pct","Bases masked (%)")]:
        vals = [fmt(rd.get(a,{}).get(key,"—")) if key in ("total_length","masked_bp","n_sequences")
                else rd.get(a,{}).get(key,"—") for a in aids]
        rr.append([label] + vals)
    lv = [", ".join([link(f"repeats/{a}.tbl","summary"),
                     link(f"repeats/{a}.buildSummary.tbl","families"),
                     link(f"repeats/{a}.repeat_landscape.html","landscape")]) for a in aids]
    rr.append(["Reports"] + lv)
    W.append(table(rh, rr))
    W.append("")
    W.append(plot_grid(aids, pd, md_pd, "Per-scaffold masking", "{aid}_repeat_masking_per_scaffold"))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 5. Non-coding RNA
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 5. Non-coding RNA\n")
    nh = ["Metric"] + [f"`{a}`" for a in aids]
    nr = [
        ["tRNA genes"]           + [fmt(met[a].get("n_trna","—")) for a in aids],
        ["ncRNA features (Rfam)"]+ [fmt(met[a].get("n_ncrna","—")) for a in aids],
        ["Reports"]              + [link("ncRNA/","ncRNA output") for _ in aids],
    ]
    W.append(table(nh, nr))
    W.append("")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("tRNA isotypes", "{aid}_trna_isotypes"),
        ("Rfam families", "{aid}_ncrna_rfam"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 6. Telomeres & Centromeres
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 6. Telomeres & Centromeres\n")
    W.append(metrics_table(aids, met, [
        ("n_telomere",   "Telomeric regions",   True),
        ("n_centromere", "Centromere features", True),
    ]))
    W.append("")
    W.append(plot_grid_multi(aids, pd, md_pd, [
        ("Telomere positions", "{aid}_telomere_positions"),
        ("Telomere density",   "{aid}_telomere_density"),
    ]))
    W.append("")

    # ═══════════════════════════════════════════════════════════════════════
    # 7. Completeness & Quality
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 7. Completeness & Quality\n")

    W.append("### BUSCO\n")
    W.append(metrics_table(aids, met, [
        ("busco_complete",   "Complete (%)",      False),
        ("busco_single",     "— Single (%)",      False),
        ("busco_duplicated", "— Duplicated (%)",  False),
        ("busco_fragmented", "Fragmented (%)",    False),
        ("busco_missing",    "Missing (%)",       False),
        ("busco_total",      "Total BUSCOs",      True),
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
    W.append(f"> Full reports: {link('qc/busco/','BUSCO output')}\n")

    W.append("\n### OMArk\n")
    oi = {}
    for aid in aids:
        of = find_file(args.omark_summaries, aid, [".summary.txt"])
        if of: oi[aid] = parse_omark_summary(read_text(of))
    if oi:
        first = next(iter(oi.values()))
        W.append(f"**Reference clade:** {first.get('clade','—')}\n")

    W.append("\n#### Completeness (conserved genes)\n")
    h = ["Metric"] + [f"`{a}`" for a in aids]
    cr = [["Conserved HOGs"] + [fmt(oi.get(a,{}).get("n_hogs","—")) for a in aids]]
    for key, label in [("omark_single","Single (%)"), ("omark_duplicated","Duplicated total (%)"),
                        ("omark_dup_unexpected","— Unexpected (%)"), ("omark_dup_expected","— Expected (%)"),
                        ("omark_missing","Missing (%)")]:
        cr.append([label] + [met[a].get(key,"—") for a in aids])
    W.append(table(h, cr))
    W.append(f"\n{img(pd, md_pd, 'C_omark_completeness.png')}\n")

    W.append("#### Consistency (proteome quality)\n")
    cr2 = [["Total proteins"] + [fmt(oi.get(a,{}).get("n_proteins","—")) for a in aids]]
    for key, label in [("omark_consistent","Consistent (%)"), ("omark_partial","— Partial hits (%)"),
                        ("omark_fragmented","— Fragmented (%)"), ("omark_inconsistent","Inconsistent (%)"),
                        ("omark_contaminant","Contaminant (%)"), ("omark_unknown","Unknown (%)")]:
        cr2.append([label] + [met[a].get(key,"—") for a in aids])
    W.append(table(h, cr2))
    W.append(f"\n{img(pd, md_pd, 'C_omark_consistency.png')}\n")
    W.append(f"> Full reports: {link('qc/omark/','OMArk output')}\n")

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
                v = ad.get(a,{}).get(label,"—")
                try:
                    vals.append(fmt(v) if "." not in v and "%" not in label.lower() else v)
                except: vals.append(v)
            ar.append([label] + vals)
        W.append(table(ah, ar))
    W.append(f"\n> Full reports: {link('qc/','QC output')}\n")

    # ═══════════════════════════════════════════════════════════════════════
    # 8. Functional Annotation
    # ═══════════════════════════════════════════════════════════════════════
    W.append("\n---\n## 8. Functional Annotation\n")
    W.append("### Tool coverage\n")
    W.append(metrics_table(aids, met, [
        ("n_proteins",       "Total proteins",    True),
        ("n_ipr",            "InterProScan hits", True),
        ("pct_ipr",          "InterProScan (%)",  False),
        ("n_eggnog",         "eggNOG hits",       True),
        ("pct_eggnog",       "eggNOG (%)",        False),
        ("n_diamond",        "DIAMOND hits",      True),
        ("pct_diamond",      "DIAMOND (%)",       False),
        ("n_kofamscan",      "KofamScan hits",    True),
        ("pct_kofamscan",    "KofamScan (%)",     False),
        ("n_deepsig",        "DeepSig hits",      True),
        ("pct_deepsig",      "DeepSig (%)",       False),
        ("n_named",          "Named products",    True),
        ("n_hypothetical",   "Hypothetical",      True),
        ("pct_named",        "Named (%)",         False),
        ("pct_hypothetical", "Hypothetical (%)",  False),
        ("unique_go_terms",  "Unique GO terms",   True),
        ("unique_kegg_kos",  "Unique KEGG KOs",   True),
    ]))
    W.append("")
    W.append("### Details\n")

    te_img = img(pd, md_pd, "C_te_flagged.png")
    if te_img:
        W.append("#### TE-derived gene models\n")
        W.append("Gene models flagged as putative TE-derived based on repeat overlap "
                 "(≥50% covered by RepeatMasker) and/or protein homology to known TE sequences. "
                 "These are retained in the GFF3 with `Note=putative_TE_derived` — "
                 "filter with `grep -v putative_TE_derived` for a clean gene set.\n")
        W.append(f"\n{te_img}\n")

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
    # 9. Tuning Recommendations
    # ═══════════════════════════════════════════════════════════════════════
    W.append("---\n## 9. Tuning Recommendations\n")
    any_flags = False
    for aid in aids:
        row = met[aid]
        flags = []
        checks = [
            ("pct_hypothetical",        50,   None, "High hypothetical rate ({v:.1f}%) — consider RNA-seq evidence or broader protein DBs."),
            ("busco_complete",          None,  80,  "Low BUSCO completeness ({v:.1f}%) — check repeat masking, contiguity, or BRAKER3 rounds."),
            ("omark_unknown",           25,   None, "High OMArk unknown rate ({v:.1f}%) — possible over-prediction. Consider TSEBRA."),
            ("omark_fragmented",        10,   None, "High OMArk fragmented rate ({v:.1f}%) — many genes shorter than expected."),
            ("omark_dup_unexpected",    15,   None, "High unexpected duplication ({v:.1f}%) — possible haplotig duplication."),
            ("repeat_masked_pct",       None,  20,  "Low repeat masking ({v:.1f}%) — fish genomes typically 30–50%."),
            ("pct_diamond",             None,  30,  "Low SwissProt coverage ({v:.1f}%) — consider TrEMBL/NR."),
            ("pct_te_flagged",           5,   None, "TE-flagged genes at {v:.1f}% — consider filtering the gene set."),
            ("mean_aed",             "Mean AED",                False),
            ("median_aed",           "Median AED",              False),
            ("pct_aed_lt05",         "AED < 0.5 (%)",           False),
            ("pct_genes_fail",           5,   None, "FAIL quality flag at {v:.1f}% of genes — inspect internal stops and splice sites."),
            ("pct_mrnas_with_internal_stop", 1, None, "Internal stop codons in {v:.1f}% of mRNAs — likely assembly or annotation errors."),
            ("pct_canonical_splice_sites", None, 96, "Low canonical splice site rate ({v:.1f}%) — check assembly quality or repeat masking completeness."),
            ("pct_single_exon_mrna",    35,   None, "High single-exon rate ({v:.1f}%) — possible TE over-prediction. Review repeat library."),
            ("pct_high_confidence",     None,  50,  "Only {v:.1f}% of genes are high-confidence — consider tightening filtering thresholds."),
            ("pct_aed_no_evidence",     20,   None, "High fraction ({v:.1f}%) of genes with no protein evidence (AED=1) — consider broader protein databases or RNA-seq."),
            ("pct_aed_lt05",            None,  50,  "Only {v:.1f}% of genes have AED < 0.5 — gene models may be poorly supported by protein evidence."),
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