/*
========================================================================================
    MERGE FUNCTIONAL ANNOTATIONS
========================================================================================
    Combine per-assembly results from InterProScan, eggNOG, DIAMOND, KofamScan,
    and DeepSig into a single unified TSV keyed on protein ID.

    Each source contributes selected columns; proteins missing from a given tool
    get "-" placeholders. The merge is a full outer join on protein ID.

    Inputs:
        - assembly_id:       Sample identifier
        - interproscan_tsv:  Merged InterProScan TSV
        - eggnog_tsv:        Merged eggNOG annotations
        - diamond_tsv:       Merged DIAMOND TSV (with header)
        - kofamscan_tsv:     Merged KofamScan TSV
        - deepsig_gff3:      Merged DeepSig GFF3

    Outputs:
        - merged_tsv:  ${assembly_id}.functional_annotations.tsv

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process MERGE_FUNCTIONAL_ANNOTATIONS {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/functional_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id),
          path(interproscan_tsv),
          path(eggnog_tsv),
          path(diamond_tsv),
          path(kofamscan_tsv),
          path(deepsig_gff3)

    output:
    tuple val(assembly_id),
          path("${assembly_id}.functional_annotations.tsv"),
          emit: merged_tsv

    script:
    """
    #!/usr/bin/env python3
    import csv, sys
    from collections import defaultdict, OrderedDict

    # ── Storage: protein_id → dict of annotation fields ──────────────────────
    annot = defaultdict(lambda: {})

    # ── 1. InterProScan (headerless TSV, 11-15 columns) ─────────────────────
    #   Col layout (0-indexed):
    #     0  protein_id
    #     3  analysis (e.g. Pfam, PANTHER)
    #     4  signature_accession
    #     5  signature_description
    #     11 IPR accession  (may be absent)
    #     12 IPR description (may be absent)
    #     13 GO terms        (may be absent)
    with open("${interproscan_tsv}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            pid = cols[0]
            rec = annot[pid]

            sig = cols[4] if len(cols) > 4 else "-"
            sig_desc = cols[5] if len(cols) > 5 else "-"
            ipr = cols[11] if len(cols) > 11 and cols[11] else "-"
            ipr_desc = cols[12] if len(cols) > 12 and cols[12] else "-"
            go = cols[13] if len(cols) > 13 and cols[13] else "-"

            # Accumulate (a protein may have many domain hits)
            entry = f"{sig}:{sig_desc}"
            rec.setdefault("ipr_signatures", []).append(entry)
            if ipr != "-":
                rec.setdefault("ipr_accessions", set()).add(ipr)
            if ipr_desc != "-":
                rec.setdefault("ipr_descriptions", set()).add(ipr_desc)
            if go != "-":
                for g in go.split("|"):
                    rec.setdefault("ipr_go", set()).add(g.strip())

    # ── 2. eggNOG (header starts with #, data lines are tab-separated) ──────
    #   Standard emapper columns (0-indexed):
    #     0  query
    #     4  eggNOG_OGs
    #     5  max_annot_lvl
    #     7  Description
    #     8  Preferred_name
    #     9  GOs
    #     11 KEGG_ko
    #     12 KEGG_Pathway
    #     15 CAZy
    #     20 PFAMs
    with open("${eggnog_tsv}") as fh:
        header_seen = False
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header_seen = True
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            pid = cols[0]
            rec = annot[pid]
            rec["eggnog_OGs"]       = cols[4]  if len(cols) > 4  and cols[4]  else "-"
            rec["eggnog_desc"]      = cols[7]  if len(cols) > 7  and cols[7]  else "-"
            rec["eggnog_name"]      = cols[8]  if len(cols) > 8  and cols[8]  else "-"
            rec["eggnog_GO"]        = cols[9]  if len(cols) > 9  and cols[9]  else "-"
            rec["eggnog_KEGG_ko"]   = cols[11] if len(cols) > 11 and cols[11] else "-"
            rec["eggnog_KEGG_path"] = cols[12] if len(cols) > 12 and cols[12] else "-"
            rec["eggnog_CAZy"]      = cols[15] if len(cols) > 15 and cols[15] else "-"
            rec["eggnog_PFAMs"]     = cols[20] if len(cols) > 20 and cols[20] else "-"

    # ── 3. DIAMOND (first line is header) ────────────────────────────────────
    #   qseqid  sseqid  stitle  pident  length  mismatch  gapopen
    #   qstart  qend  sstart  send  evalue  bitscore
    with open("${diamond_tsv}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            pid = row.get("qseqid", "").strip()
            if not pid:
                continue
            rec = annot[pid]
            # Keep best hit (first encountered per protein since DIAMOND is pre-sorted)
            if "diamond_hit" not in rec:
                rec["diamond_hit"]    = row.get("sseqid", "-")
                rec["diamond_desc"]   = row.get("stitle", "-")
                rec["diamond_pident"] = row.get("pident", "-")
                rec["diamond_evalue"] = row.get("evalue", "-")

    # ── 4. KofamScan (first line is header, detail-tsv format) ───────────────
    #   Fields: significance, gene_name, KO, thrshld, score, E-value, KO_definition
    #   '*' in significance column means above threshold
    with open("${kofamscan_tsv}") as fh:
        for i, line in enumerate(fh):
            if i == 0 or line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            # detail-tsv: col0=significance marker, col1=gene_name, col2=KO, ...
            if len(cols) < 7:
                continue
            sig_marker = cols[0].strip()
            pid = cols[1].strip()
            ko = cols[2].strip() if cols[2].strip() else "-"
            score = cols[4].strip() if len(cols) > 4 else "-"
            evalue = cols[5].strip() if len(cols) > 5 else "-"
            ko_def = cols[6].strip() if len(cols) > 6 else "-"

            rec = annot[pid]
            # Keep only significant hits (marked with '*'), or best hit
            if "kofam_KO" not in rec or sig_marker == "*":
                rec["kofam_KO"]     = ko
                rec["kofam_score"]  = score
                rec["kofam_evalue"] = evalue
                rec["kofam_def"]    = ko_def
                rec["kofam_sig"]    = sig_marker if sig_marker else "-"

    # ── 5. DeepSig (GFF3 format) ────────────────────────────────────────────
    #   Attributes field contains signal peptide prediction info
    with open("${deepsig_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            pid = cols[0]
            feat_type = cols[2] if len(cols) > 2 else "-"
            start = cols[3] if len(cols) > 3 else "-"
            end = cols[4] if len(cols) > 4 else "-"
            score = cols[5] if len(cols) > 5 else "-"
            attributes = cols[8] if len(cols) > 8 else "-"

            rec = annot[pid]
            if "deepsig_type" not in rec:
                rec["deepsig_type"]  = feat_type
                rec["deepsig_start"] = start
                rec["deepsig_end"]   = end
                rec["deepsig_score"] = score
                rec["deepsig_attr"]  = attributes

    # ── Build output ─────────────────────────────────────────────────────────
    out_columns = [
        "protein_id",
        # InterProScan
        "ipr_signatures", "ipr_accessions", "ipr_descriptions", "ipr_go",
        # eggNOG
        "eggnog_OGs", "eggnog_name", "eggnog_desc", "eggnog_GO",
        "eggnog_KEGG_ko", "eggnog_KEGG_path", "eggnog_CAZy", "eggnog_PFAMs",
        # DIAMOND
        "diamond_hit", "diamond_desc", "diamond_pident", "diamond_evalue",
        # KofamScan
        "kofam_KO", "kofam_def", "kofam_score", "kofam_evalue", "kofam_sig",
        # DeepSig
        "deepsig_type", "deepsig_start", "deepsig_end", "deepsig_score", "deepsig_attr",
    ]

    with open("${assembly_id}.functional_annotations.tsv", "w") as out:
        out.write("\\t".join(out_columns) + "\\n")

        for pid in sorted(annot.keys()):
            rec = annot[pid]
            vals = [pid]

            # InterProScan — collapse accumulated lists/sets
            vals.append("; ".join(rec.get("ipr_signatures", ["-"])))
            vals.append("; ".join(sorted(rec.get("ipr_accessions", {"-"}))))
            vals.append("; ".join(sorted(rec.get("ipr_descriptions", {"-"}))))
            vals.append("; ".join(sorted(rec.get("ipr_go", {"-"}))))

            # Scalar fields
            for col in out_columns[5:]:
                vals.append(rec.get(col, "-"))

            out.write("\\t".join(str(v) for v in vals) + "\\n")

    # Quick stats
    total = len(annot)
    has_ipr = sum(1 for r in annot.values() if "ipr_signatures" in r)
    has_egg = sum(1 for r in annot.values() if "eggnog_OGs" in r)
    has_dia = sum(1 for r in annot.values() if "diamond_hit" in r)
    has_kof = sum(1 for r in annot.values() if "kofam_KO" in r)
    has_sig = sum(1 for r in annot.values() if "deepsig_type" in r)
    print(f"Merged {total} proteins: IPR={has_ipr} eggNOG={has_egg} DIAMOND={has_dia} KofamScan={has_kof} DeepSig={has_sig}", file=sys.stderr)
    """

    stub:
    """
    touch ${assembly_id}.functional_annotations.tsv
    """
}