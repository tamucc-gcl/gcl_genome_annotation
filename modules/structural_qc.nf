/*
========================================================================================
    STRUCTURAL QC MODULE
========================================================================================
    Assembly-level structural annotation quality assessment.

    Reads the enriched functional annotations TSV (which now contains per-protein
    quality metrics including AED scores) and the TE-flagged BRAKER GFF3 to produce:

      1. Per-gene quality summary (worst-transcript flag, TE status, exon counts,
         AED, etc.)
      2. Assembly-level QC statistics table (includes AED distribution stats)
      3. High-confidence gene list (PASS, non-TE-flagged)

    The enriched TSV uses original BRAKER protein/mRNA IDs (e.g. g1.t1).
    The flagged GFF3 uses the same IDs. The high-confidence gene list is therefore
    in original BRAKER gene IDs (e.g. g1), which correspond 1:1 with the IDs
    in ANNOTATE_GFF3's renamed output (${assembly_id}_g1).

    Inputs:
        - assembly_id:     Sample identifier
        - functional_tsv:  Enriched functional annotations TSV
                           (from MERGE_FUNCTIONAL_ANNOTATIONS; includes aed column)
        - flagged_gff3:    TE-flagged BRAKER GFF3 (from FLAG_TE_GENES)

    Outputs:
        - structural_qc:       ${assembly_id}.structural_qc.tsv        (assembly-level stats)
        - structural_qc_genes: ${assembly_id}.structural_qc_genes.tsv  (per-gene summary)
        - high_confidence:     ${assembly_id}.high_confidence_genes.txt (PASS gene IDs)

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process STRUCTURAL_QC {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/qc/structural", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(functional_tsv), path(flagged_gff3)

    output:
    tuple val(assembly_id),
          path("${assembly_id}.structural_qc.tsv"),
          emit: structural_qc
    tuple val(assembly_id),
          path("${assembly_id}.structural_qc_genes.tsv"),
          emit: structural_qc_genes
    tuple val(assembly_id),
          path("${assembly_id}.high_confidence_genes.txt"),
          emit: high_confidence

    script:
    """
    #!/usr/bin/env python3
    import csv, sys, statistics
    from collections import defaultdict

    ASM = "${assembly_id}"

    # ── 1. Parse GFF3: build gene → mRNA hierarchy and TE-flag status ────────
    gene_to_mrnas = defaultdict(list)   # gene_id → [mrna_id, ...]
    mrna_to_gene  = {}                   # mrna_id → gene_id
    te_flagged_genes = set()             # gene IDs with putative_TE_derived note

    with open("${flagged_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            ftype  = cols[2]
            attrs  = {}
            for tok in cols[8].split(";"):
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attrs[k.strip()] = v.strip()

            fid    = attrs.get("ID", "")
            parent = attrs.get("Parent", "")

            if ftype == "gene":
                if "putative_TE_derived" in cols[8]:
                    te_flagged_genes.add(fid)

            elif ftype in ("mRNA", "transcript") and parent:
                gene_to_mrnas[parent].append(fid)
                mrna_to_gene[fid] = parent

    # ── 2. Load functional TSV: per-mRNA quality metrics ─────────────────────
    FLAG_ORDER = {"FAIL": 2, "WARN": 1, "PASS": 0}

    mrna_metrics = {}   # mrna_id → row dict

    with open("${functional_tsv}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            pid = row.get("protein_id", "").strip()
            if pid:
                mrna_metrics[pid] = row

    # ── 3. Compute per-gene summary ───────────────────────────────────────────

    def safe_int(val, default=None):
        try: return int(val)
        except (ValueError, TypeError): return default

    def safe_float(val, default=None):
        try: return float(val)
        except (ValueError, TypeError): return default

    def worst_flag(*flags):
        best_rank = -1
        worst = "PASS"
        for f in flags:
            if f in FLAG_ORDER and FLAG_ORDER[f] > best_rank:
                best_rank = FLAG_ORDER[f]
                worst = f
        return worst

    gene_rows = []

    all_gene_ids = set(gene_to_mrnas.keys())
    for mid, gid in mrna_to_gene.items():
        all_gene_ids.add(gid)

    for gene_id in sorted(all_gene_ids):
        mrna_ids = gene_to_mrnas.get(gene_id, [])

        flags          = []
        prot_lens      = []
        exon_counts    = []
        intron_counts  = []
        pct_can_vals   = []
        n_short_all    = 0
        n_noncan_all   = 0
        has_func_any   = False
        omark_classes  = []
        aed_vals_gene  = []
        reasons_all    = []

        for mid in mrna_ids:
            m = mrna_metrics.get(mid, {})
            flags.append(m.get("quality_flag", "PASS"))
            reasons = m.get("quality_reasons", ".")
            if reasons and reasons != ".":
                reasons_all.append(f"{mid}:{reasons}")

            plen = safe_int(m.get("protein_length_aa"))
            if plen is not None:
                prot_lens.append(plen)

            ne = safe_int(m.get("n_exons"))
            if ne is not None:
                exon_counts.append(ne)

            ni = safe_int(m.get("n_introns"))
            if ni is not None:
                intron_counts.append(ni)

            pct_can = safe_float(m.get("pct_canonical_splice"))
            if pct_can is not None:
                pct_can_vals.append(pct_can)

            ns = safe_int(m.get("n_short_introns_lt60"), 0)
            n_short_all += ns

            nnc = safe_int(m.get("n_noncanonical_introns"), 0)
            n_noncan_all += nnc

            has_func = any(
                m.get(k, "-") not in ("-", "", None)
                for k in ("ipr_signatures", "eggnog_OGs", "diamond_hit", "kofam_KO")
            )
            if has_func:
                has_func_any = True

            ok = m.get("omark_class", "-")
            if ok and ok not in ("-", "."):
                omark_classes.append(ok)

            aed_val = safe_float(m.get("aed"))
            if aed_val is not None:
                aed_vals_gene.append(aed_val)

        n_mrnas  = len(mrna_ids)
        n_pass   = flags.count("PASS")
        n_warn   = flags.count("WARN")
        n_fail   = flags.count("FAIL")
        wflag    = worst_flag(*flags) if flags else "PASS"
        bflag    = min(flags, key=lambda f: FLAG_ORDER.get(f, 0)) if flags else "PASS"
        is_te    = "Y" if gene_id in te_flagged_genes else "N"

        mean_plen = f"{statistics.mean(prot_lens):.1f}" if prot_lens else "-"
        mean_exon = f"{statistics.mean(exon_counts):.1f}" if exon_counts else "-"
        min_pct_can = f"{min(pct_can_vals):.1f}" if pct_can_vals else "."
        best_omark = (
            "Consistent" if "Consistent" in omark_classes else
            omark_classes[0] if omark_classes else "-"
        )

        min_aed = f"{min(aed_vals_gene):.4f}" if aed_vals_gene else "."
        mean_aed_gene = f"{statistics.mean(aed_vals_gene):.4f}" if aed_vals_gene else "."

        # HC criteria: worst transcript flag is PASS, not TE-flagged
        is_hc = (wflag == "PASS" and is_te == "N")

        gene_rows.append({
            "gene_id":            gene_id,
            "n_mrnas":            n_mrnas,
            "worst_quality_flag": wflag,
            "best_quality_flag":  bflag,
            "n_pass":             n_pass,
            "n_warn":             n_warn,
            "n_fail":             n_fail,
            "is_te_flagged":      is_te,
            "is_high_confidence": "Y" if is_hc else "N",
            "mean_protein_len_aa":mean_plen,
            "mean_exons_per_mrna":mean_exon,
            "min_pct_canonical_splice": min_pct_can,
            "n_short_introns_lt60":     n_short_all,
            "n_noncanonical_introns":   n_noncan_all,
            "has_functional_support":   "Y" if has_func_any else "N",
            "best_omark_class":         best_omark,
            "min_aed":                  min_aed,
            "mean_aed":                 mean_aed_gene,
            "quality_reasons":          "|".join(reasons_all) if reasons_all else ".",
        })

    # ── 4. Write per-gene TSV ─────────────────────────────────────────────────
    gene_cols = [
        "gene_id", "n_mrnas",
        "worst_quality_flag", "best_quality_flag",
        "n_pass", "n_warn", "n_fail",
        "is_te_flagged", "is_high_confidence",
        "mean_protein_len_aa", "mean_exons_per_mrna",
        "min_pct_canonical_splice", "n_short_introns_lt60",
        "n_noncanonical_introns", "has_functional_support",
        "best_omark_class", "min_aed", "mean_aed",
        "quality_reasons",
    ]

    with open(f"{ASM}.structural_qc_genes.tsv", "w") as out:
        out.write("\\t".join(gene_cols) + "\\n")
        for row in gene_rows:
            out.write("\\t".join(str(row.get(c, "-")) for c in gene_cols) + "\\n")

    # ── 5. Write high-confidence gene list ────────────────────────────────────
    hc_genes = [r["gene_id"] for r in gene_rows if r["is_high_confidence"] == "Y"]
    with open(f"{ASM}.high_confidence_genes.txt", "w") as out:
        out.write("# High-confidence genes: PASS quality flag, not TE-flagged\\n")
        out.write(f"# Assembly: {ASM}\\n")
        out.write(f"# Total: {len(hc_genes)}\\n")
        for gid in hc_genes:
            out.write(gid + "\\n")

    # ── 6. Compute assembly-level QC statistics ───────────────────────────────
    n_genes  = len(gene_rows)
    n_mrnas_all = sum(r["n_mrnas"] for r in gene_rows)
    n_g_pass = sum(1 for r in gene_rows if r["worst_quality_flag"] == "PASS")
    n_g_warn = sum(1 for r in gene_rows if r["worst_quality_flag"] == "WARN")
    n_g_fail = sum(1 for r in gene_rows if r["worst_quality_flag"] == "FAIL")
    n_te     = sum(1 for r in gene_rows if r["is_te_flagged"] == "Y")
    n_hc     = len(hc_genes)
    n_func   = sum(1 for r in gene_rows if r["has_functional_support"] == "Y")

    # Per-mRNA stats from enriched TSV
    prot_lens_all, exon_counts_all, intron_lens_all = [], [], []
    pct_can_all, n_single_exon = [], 0
    n_start_m, n_internal_stop, n_no_stop = 0, 0, 0
    n_total_introns, n_total_noncan, n_total_short = 0, 0, 0
    omark_dist = defaultdict(int)

    for m in mrna_metrics.values():
        plen = safe_int(m.get("protein_length_aa"))
        if plen is not None:
            prot_lens_all.append(plen)

        ne = safe_int(m.get("n_exons"))
        if ne is not None:
            exon_counts_all.append(ne)
            if ne == 1:
                n_single_exon += 1

        ni = safe_int(m.get("n_introns"), 0)
        n_total_introns += ni

        nnc = safe_int(m.get("n_noncanonical_introns"), 0)
        n_total_noncan += nnc

        ns = safe_int(m.get("n_short_introns_lt60"), 0)
        n_total_short += ns

        pct_can = safe_float(m.get("pct_canonical_splice"))
        if pct_can is not None:
            pct_can_all.append(pct_can)

        if m.get("has_start_M") == "Y":
            n_start_m += 1
        if m.get("has_internal_stop") == "Y":
            n_internal_stop += 1
        if m.get("has_stop_codon") == "N":
            n_no_stop += 1

        ok = m.get("omark_class", "-")
        if ok and ok not in ("-", "."):
            omark_dist[ok] += 1

    n_mrnas_tot = len(mrna_metrics)
    pct = lambda n, d: f"{100*n/d:.2f}" if d > 0 else "0.00"

    def _med(lst):
        return f"{statistics.median(lst):.1f}" if lst else "-"
    def _mean(lst):
        return f"{statistics.mean(lst):.1f}" if lst else "-"

    omark_str = ",".join(f"{k}:{v}" for k, v in sorted(omark_dist.items()))

    # ── AED distribution (from enriched functional TSV) ───────────────────────
    #   The aed column may be "." for all rows if no protein DB was available;
    #   in that case all_aed will be empty and stats will be ".".
    all_aed = []
    for m in mrna_metrics.values():
        a = safe_float(m.get("aed"))
        if a is not None:
            all_aed.append(a)

    n_aed_lt03 = sum(1 for a in all_aed if a < 0.3)
    n_aed_lt05 = sum(1 for a in all_aed if a < 0.5)
    n_aed_eq1  = sum(1 for a in all_aed if a >= 0.999)

    # ── Build stats dict ──────────────────────────────────────────────────────
    stats = {
        "assembly_id":                ASM,
        "n_genes":                    n_genes,
        "n_mrnas":                    n_mrnas_tot,
        "n_genes_pass":               n_g_pass,
        "n_genes_warn":               n_g_warn,
        "n_genes_fail":               n_g_fail,
        "pct_genes_pass":             pct(n_g_pass, n_genes),
        "pct_genes_warn":             pct(n_g_warn, n_genes),
        "pct_genes_fail":             pct(n_g_fail, n_genes),
        "n_te_flagged":               n_te,
        "pct_te_flagged":             pct(n_te, n_genes),
        "n_high_confidence":          n_hc,
        "pct_high_confidence":        pct(n_hc, n_genes),
        "n_with_functional_support":  n_func,
        "pct_with_functional_support":pct(n_func, n_genes),
        "n_single_exon_mrna":         n_single_exon,
        "pct_single_exon_mrna":       pct(n_single_exon, n_mrnas_tot),
        "n_mrnas_with_start_m":       n_start_m,
        "pct_mrnas_with_start_m":     pct(n_start_m, n_mrnas_tot),
        "n_mrnas_with_internal_stop": n_internal_stop,
        "pct_mrnas_with_internal_stop":pct(n_internal_stop, n_mrnas_tot),
        "n_mrnas_no_stop_codon":      n_no_stop,
        "pct_mrnas_no_stop_codon":    pct(n_no_stop, n_mrnas_tot),
        "mean_protein_length_aa":     _mean(prot_lens_all),
        "median_protein_length_aa":   _med(prot_lens_all),
        "mean_exons_per_mrna":        _mean(exon_counts_all),
        "median_exons_per_mrna":      _med(exon_counts_all),
        "n_total_introns":            n_total_introns,
        "n_total_noncanonical":       n_total_noncan,
        "pct_canonical_splice_sites": pct(n_total_introns - n_total_noncan, n_total_introns),
        "n_short_introns_lt60":       n_total_short,
        "pct_short_introns_lt60":     pct(n_total_short, n_total_introns),
        "mean_pct_canonical_per_mrna":_mean(pct_can_all),
        "omark_class_distribution":   omark_str if omark_str else "-",
        # AED stats — will be "." when no protein DB was available
        "mean_aed":                   f"{statistics.mean(all_aed):.4f}" if all_aed else ".",
        "median_aed":                 f"{statistics.median(all_aed):.4f}" if all_aed else ".",
        "pct_aed_lt03":               pct(n_aed_lt03, len(all_aed)) if all_aed else ".",
        "pct_aed_lt05":               pct(n_aed_lt05, len(all_aed)) if all_aed else ".",
        "pct_aed_no_evidence":        pct(n_aed_eq1, len(all_aed)) if all_aed else ".",
    }

    stat_cols = list(stats.keys())
    with open(f"{ASM}.structural_qc.tsv", "w") as out:
        out.write("\\t".join(stat_cols) + "\\n")
        out.write("\\t".join(str(stats[c]) for c in stat_cols) + "\\n")

    # ── 7. Summary to stderr ─────────────────────────────────────────────────
    print(f"Structural QC [{ASM}]:", file=sys.stderr)
    print(f"  Genes: {n_genes} total | PASS={n_g_pass} WARN={n_g_warn} FAIL={n_g_fail}", file=sys.stderr)
    print(f"  High-confidence: {n_hc} ({pct(n_hc, n_genes)}%)", file=sys.stderr)
    print(f"  Splice sites canonical: {pct(n_total_introns - n_total_noncan, n_total_introns)}%", file=sys.stderr)
    print(f"  Single-exon mRNAs: {pct(n_single_exon, n_mrnas_tot)}%", file=sys.stderr)
    print(f"  Internal stops: {pct(n_internal_stop, n_mrnas_tot)}%", file=sys.stderr)
    if all_aed:
        print(f"  AED: mean={statistics.mean(all_aed):.4f}, "
              f"<0.3={pct(n_aed_lt03, len(all_aed))}%, "
              f"<0.5={pct(n_aed_lt05, len(all_aed))}%, "
              f"no_ev={pct(n_aed_eq1, len(all_aed))}%", file=sys.stderr)
    else:
        print(f"  AED: not computed (no protein DB)", file=sys.stderr)
    """

    stub:
    """
    touch ${assembly_id}.structural_qc.tsv
    touch ${assembly_id}.structural_qc_genes.tsv
    touch ${assembly_id}.high_confidence_genes.txt
    """
}