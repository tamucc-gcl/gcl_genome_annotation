/*
========================================================================================
    MERGE FUNCTIONAL ANNOTATIONS
========================================================================================
    Combine per-assembly results from InterProScan, eggNOG, DIAMOND, KofamScan,
    DeepSig, and OMArk into a single unified TSV keyed on protein ID.

    Also computes per-protein quality metrics:
      - Protein sequence integrity (start M, internal stops, length)
      - Structural metrics from GFF3 (exon/intron counts, CDS length)
      - Splice-site canonicality from genome sequence (GT-AG audit)
      - OMArk per-protein classification (Consistent/Partial/etc.)
      - Composite quality flag (PASS / WARN / FAIL)

    Inputs:
        - assembly_id:       Sample identifier
        - interproscan_tsv:  Merged InterProScan TSV
        - eggnog_tsv:        Merged eggNOG annotations
        - diamond_tsv:       Merged DIAMOND TSV (with header)
        - kofamscan_tsv:     Merged KofamScan TSV
        - deepsig_gff3:      Merged DeepSig GFF3
        - proteins_fasta:    BRAKER protein FASTA (for sequence integrity)
        - braker_gff3:       TE-flagged BRAKER GFF3 (for structural metrics)
        - genome_fasta:      Filtered genome assembly (for splice-site lookup)
        - omark_ump:         OMArk Unified Model Predictions (per-protein classification)

    Outputs:
        - merged_tsv:  ${assembly_id}.functional_annotations.tsv

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process MERGE_FUNCTIONAL_ANNOTATIONS {
    tag "${assembly_id}"
    label 'merge_functional_annotations'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id),
          path(interproscan_tsv),
          path(eggnog_tsv),
          path(diamond_tsv),
          path(kofamscan_tsv),
          path(deepsig_gff3),
          path(proteins_fasta),
          path(braker_gff3),
          path(genome_fasta),
          path(omark_ump),
          path(aed_scores) 

    output:
    tuple val(assembly_id),
          path("${assembly_id}.functional_annotations.tsv"),
          emit: merged_tsv

    script:
    """
    #!/usr/bin/env python3
    import csv, re, sys, statistics
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
    with open("${eggnog_tsv}") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
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
    with open("${diamond_tsv}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            pid = row.get("qseqid", "").strip()
            if not pid:
                continue
            rec = annot[pid]
            if "diamond_hit" not in rec:
                rec["diamond_hit"]    = row.get("sseqid", "-")
                rec["diamond_desc"]   = row.get("stitle", "-")
                rec["diamond_pident"] = row.get("pident", "-")
                rec["diamond_evalue"] = row.get("evalue", "-")

    # ── 4. KofamScan (first line is header, detail-tsv format) ───────────────
    with open("${kofamscan_tsv}") as fh:
        for i, line in enumerate(fh):
            if i == 0 or line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 7:
                continue
            sig_marker = cols[0].strip()
            pid = cols[1].strip()
            ko = cols[2].strip() if cols[2].strip() else "-"
            score = cols[4].strip() if len(cols) > 4 else "-"
            evalue = cols[5].strip() if len(cols) > 5 else "-"
            ko_def = cols[6].strip() if len(cols) > 6 else "-"

            rec = annot[pid]
            if "kofam_KO" not in rec or sig_marker == "*":
                rec["kofam_KO"]     = ko
                rec["kofam_score"]  = score
                rec["kofam_evalue"] = evalue
                rec["kofam_def"]    = ko_def
                rec["kofam_sig"]    = sig_marker if sig_marker else "-"

    # ── 5. DeepSig (GFF3 format) ────────────────────────────────────────────
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

    # ── 6. Pre-populate annot from GFF3 to ensure all mRNAs are represented ──
    #   Proteins with no functional tool hits still need structural metrics.
    with open("${braker_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9 or cols[2] not in ("mRNA", "transcript"):
                continue
            attrs = {}
            for tok in cols[8].split(";"):
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attrs[k.strip()] = v.strip()
            fid = attrs.get("ID", "")
            if fid:
                _ = annot[fid]   # ensure entry exists

    # ── 7. Protein sequence metrics ──────────────────────────────────────────
    #   For each protein: length, start Met, internal stops, terminal stop.
    #   BRAKER protein FASTA headers: >g1.t1 [optional description]
    def _process_prot(pid, seq):
        if not seq or not pid:
            return
        rec = annot[pid]
        rec["protein_length_aa"] = len(seq.rstrip("*"))
        rec["has_start_M"] = "Y" if seq and seq[0] == "M" else "N"
        body = seq[1:-1] if len(seq) > 2 else seq[1:]  # skip first residue and terminal *
        rec["has_internal_stop"] = "Y" if "*" in body else "N"
        rec["has_stop_codon"] = "Y" if seq.rstrip() and seq.rstrip()[-1] == "*" else "N"

    current_pid = None
    seq_parts = []
    with open("${proteins_fasta}") as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if current_pid is not None:
                    _process_prot(current_pid, "".join(seq_parts))
                current_pid = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_pid is not None:
        _process_prot(current_pid, "".join(seq_parts))

    # ── 8. Structural metrics + splice-site audit ─────────────────────────────
    #   Load genome into memory for splice-site dinucleotide lookup.
    #   Introns are inferred from exon coordinate gaps (more robust than
    #   relying on explicit intron features in the GFF3).

    # Load genome
    genome = {}
    current_scaff = None
    scaff_parts = []
    with open("${genome_fasta}") as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if current_scaff is not None:
                    genome[current_scaff] = "".join(scaff_parts).upper()
                current_scaff = line[1:].split()[0]
                scaff_parts = []
            else:
                scaff_parts.append(line)
    if current_scaff is not None:
        genome[current_scaff] = "".join(scaff_parts).upper()

    COMP = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")
    def revcomp(s):
        return s.translate(COMP)[::-1]

    def splice_class(donor, acceptor):
        d, a = donor.upper(), acceptor.upper()
        if d == "GT" and a == "AG": return "GTAG"
        if d == "GC" and a == "AG": return "GCAG"
        if d == "AT" and a == "AC": return "ATAC"
        return "non_canonical"

    # Parse GFF3: build mRNA → [(seqid, exon_start, exon_end, strand)] and CDS lengths
    mrna_meta = {}    # mrna_id → (seqid, strand)
    mrna_exons = {}   # mrna_id → [(start, end)]
    mrna_cds_bp = {}  # mrna_id → total CDS bp

    with open("${braker_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            seqid  = cols[0]
            ftype  = cols[2]
            start  = int(cols[3])
            end    = int(cols[4])
            strand = cols[6]
            attrs  = {}
            for tok in cols[8].split(";"):
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attrs[k.strip()] = v.strip()

            parent = attrs.get("Parent", "")
            fid    = attrs.get("ID", "")

            if ftype in ("mRNA", "transcript"):
                mrna_meta[fid] = (seqid, strand)
                mrna_exons.setdefault(fid, [])

            elif ftype == "exon" and parent:
                mrna_exons.setdefault(parent, []).append((start, end))

            elif ftype == "CDS" and parent:
                mrna_cds_bp[parent] = mrna_cds_bp.get(parent, 0) + (end - start + 1)

    # Compute structural metrics per mRNA
    for mrna_id, exons in mrna_exons.items():
        if not exons:
            continue
        exons_s = sorted(exons, key=lambda x: x[0])
        seqid, strand = mrna_meta.get(mrna_id, ("", "+"))
        n_exons = len(exons_s)

        rec = annot[mrna_id]
        rec["n_exons"]       = n_exons
        rec["is_single_exon"]= "Y" if n_exons == 1 else "N"
        rec["min_exon_len"]  = min(e[1] - e[0] + 1 for e in exons_s)
        rec["cds_length_bp"] = mrna_cds_bp.get(mrna_id, 0)

        if n_exons < 2:
            rec["n_introns"]               = 0
            rec["n_canonical_introns"]     = 0
            rec["n_noncanonical_introns"]  = 0
            rec["pct_canonical_splice"]    = "."
            rec["min_intron_len"]          = "."
            rec["n_short_introns_lt60"]    = 0
            rec["donor_acceptor_summary"]  = "."
            continue

        # Infer introns from exon gaps
        introns = []
        for i in range(len(exons_s) - 1):
            istart = exons_s[i][1] + 1       # 1-based intron start
            iend   = exons_s[i+1][0] - 1     # 1-based intron end
            if iend >= istart:
                introns.append((istart, iend))

        if not introns:
            rec["n_introns"] = 0
            rec["n_canonical_introns"] = 0
            rec["n_noncanonical_introns"] = 0
            rec["pct_canonical_splice"] = "."
            rec["min_intron_len"] = "."
            rec["n_short_introns_lt60"] = 0
            rec["donor_acceptor_summary"] = "."
            continue

        intron_lens = [e[1] - e[0] + 1 for e in introns]
        rec["n_introns"]            = len(introns)
        rec["min_intron_len"]       = min(intron_lens)
        rec["n_short_introns_lt60"] = sum(1 for l in intron_lens if l < 60)

        # Splice-site dinucleotide lookup
        scaff_seq = genome.get(seqid, "")
        n_can = 0
        n_noncan = 0
        class_counts = defaultdict(int)

        for istart, iend in introns:
            if not scaff_seq or iend > len(scaff_seq) or istart < 1 or iend < istart:
                n_noncan += 1
                class_counts["non_canonical"] += 1
                continue
            # Out-of-bounds guard for 2-bp lookups
            if istart + 1 > len(scaff_seq) or iend < 2:
                n_noncan += 1
                class_counts["non_canonical"] += 1
                continue

            # Donor/acceptor in transcript orientation
            if strand == "+":
                donor    = scaff_seq[istart - 1 : istart + 1]   # first 2 bp
                acceptor = scaff_seq[iend - 2   : iend]          # last 2 bp
            else:
                donor    = revcomp(scaff_seq[iend - 2  : iend])
                acceptor = revcomp(scaff_seq[istart - 1: istart + 1])

            sc = splice_class(donor, acceptor)
            class_counts[sc] += 1
            if sc != "non_canonical":
                n_can += 1
            else:
                n_noncan += 1

        total = n_can + n_noncan
        rec["n_canonical_introns"]    = n_can
        rec["n_noncanonical_introns"] = n_noncan
        rec["pct_canonical_splice"]   = f"{100 * n_can / total:.1f}" if total > 0 else "."
        # Summary: GTAG:N,GCAG:N,ATAC:N,non_canonical:N
        rec["donor_acceptor_summary"] = ",".join(
            f"{k}:{v}" for k, v in sorted(class_counts.items()) if v > 0
        ) or "."

    # ── 9. OMArk per-protein classification (.ump file) ──────────────────────
    # Format: FASTA-style, proteins grouped under >Category_Subtype headers
    # Top-level class: Consistent, Inconsistent, Contamination, Unknown
    # Subtype: Full, Partial, Fragment (absent for Unknown)
    current_class = None
    with open("${omark_ump}") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_class = line[1:]  # e.g. "Consistent_Partial"
            elif current_class is not None:
                rec = annot[line]
                rec["omark_class"] = current_class

    # ── 9b. AED scores ──────────────────────────────────────────────────────
    #   Read per-transcript AED from COMPUTE_AED output.
    #   File may be header-only if no protein DB was available.
    import csv as csv_mod
    with open("${aed_scores}") as fh:
        reader = csv_mod.DictReader(fh, delimiter="\\t")
        for row in reader:
            tid = row.get("transcript_id", "").strip()
            if tid and tid in annot:
                annot[tid]["aed"]                  = row.get("aed", ".")
                annot[tid]["aed_best_hit"]         = row.get("best_protein_hit", ".")
                annot[tid]["aed_protein_identity"] = row.get("protein_identity", ".")
                annot[tid]["aed_n_evidence_hits"]  = row.get("n_evidence_hits", ".")

    # ── 10. Quality flags ────────────────────────────────────────────────────
    #   FAIL: internal stop codon | very short protein (<30 aa) |
    #         >50% of introns are non-canonical (min 2 introns)
    #   WARN: no start Met | short protein (<100 aa) | any intron <60 bp |
    #         single-exon with no functional support |
    #         OMArk Fragmented or Inconsistent | >10% non-canonical splice
    #   PASS: all other cases

    FAIL_REASONS = []  # reused below

    for pid, rec in annot.items():
        reasons = []
        flag = "PASS"

        # ── FAIL conditions ────────────────────────────────────────────────
        if rec.get("has_internal_stop") == "Y":
            reasons.append("internal_stop_codon")
            flag = "FAIL"

        prot_len = rec.get("protein_length_aa")
        if prot_len is not None:
            try:
                prot_len = int(prot_len)
                if prot_len < 30:
                    reasons.append(f"very_short_protein:{prot_len}aa")
                    flag = "FAIL"
            except (ValueError, TypeError):
                pass

        n_introns = rec.get("n_introns", 0)
        n_noncan  = rec.get("n_noncanonical_introns", 0)
        try:
            n_introns = int(n_introns)
            n_noncan  = int(n_noncan)
            if n_introns >= 2 and n_noncan / n_introns > 0.5:
                reasons.append(f"majority_noncanonical_splice:{n_noncan}/{n_introns}")
                flag = "FAIL"
        except (ValueError, TypeError, ZeroDivisionError):
            pass

        # ── WARN conditions (only if not already FAIL) ────────────────────
        if flag != "FAIL":
            if rec.get("has_start_M") == "N":
                reasons.append("no_start_met")
                flag = "WARN"

            if prot_len is not None:
                try:
                    if int(prot_len) < 100:
                        reasons.append(f"short_protein:{prot_len}aa")
                        flag = "WARN"
                except (ValueError, TypeError):
                    pass

            n_short = rec.get("n_short_introns_lt60", 0)
            try:
                if int(n_short) > 0:
                    reasons.append(f"short_introns_lt60:{n_short}")
                    flag = "WARN"
            except (ValueError, TypeError):
                pass

            is_single = rec.get("is_single_exon") == "Y"
            has_func  = any(
                rec.get(k, "-") not in ("-", "", None)
                for k in ("ipr_signatures", "eggnog_OGs", "diamond_hit", "kofam_KO")
            )
            if is_single and not has_func:
                reasons.append("single_exon_no_functional_support")
                flag = "WARN"

            omark = rec.get("omark_class", "-")
            if omark is not None and (
                omark.startswith("Inconsistent") or
                omark.startswith("Contamination") or
                omark == "Consistent_Fragment" or
                omark == "Consistent_Partial"
            ):

            pct_can_str = rec.get("pct_canonical_splice", ".")
            if pct_can_str not in (".", "-", None, ""):
                try:
                    pct_can = float(pct_can_str)
                    if pct_can < 90.0 and n_introns >= 2:
                        reasons.append(f"low_canonical_splice_pct:{pct_can:.1f}")
                        flag = "WARN"
                except (ValueError, TypeError):
                    pass

        rec["quality_flag"]    = flag
        rec["quality_reasons"] = ";".join(reasons) if reasons else "."

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
        # Protein sequence integrity
        "protein_length_aa", "has_start_M", "has_internal_stop", "has_stop_codon",
        # Structural metrics
        "n_exons", "n_introns", "cds_length_bp", "is_single_exon",
        "min_exon_len", "min_intron_len", "n_short_introns_lt60",
        # Splice-site audit
        "n_canonical_introns", "n_noncanonical_introns",
        "pct_canonical_splice", "donor_acceptor_summary",
        # OMArk
        "omark_class",
        # AED
        "aed", "aed_best_hit", "aed_protein_identity", "aed_n_evidence_hits",
        # Quality
        "quality_flag", "quality_reasons",
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

            # All scalar fields
            for col in out_columns[5:]:
                vals.append(str(rec.get(col, "-")))

            out.write("\\t".join(str(v) for v in vals) + "\\n")

    # ── Stats ─────────────────────────────────────────────────────────────────
    total  = len(annot)
    n_ipr  = sum(1 for r in annot.values() if "ipr_signatures" in r)
    n_egg  = sum(1 for r in annot.values() if "eggnog_OGs" in r)
    n_dia  = sum(1 for r in annot.values() if "diamond_hit" in r)
    n_kof  = sum(1 for r in annot.values() if "kofam_KO" in r)
    n_sig  = sum(1 for r in annot.values() if "deepsig_type" in r)
    n_pass = sum(1 for r in annot.values() if r.get("quality_flag") == "PASS")
    n_warn = sum(1 for r in annot.values() if r.get("quality_flag") == "WARN")
    n_fail = sum(1 for r in annot.values() if r.get("quality_flag") == "FAIL")
    print(f"Merged {total} proteins: IPR={n_ipr} eggNOG={n_egg} "
          f"DIAMOND={n_dia} KofamScan={n_kof} DeepSig={n_sig}", file=sys.stderr)
    print(f"Quality flags: PASS={n_pass} WARN={n_warn} FAIL={n_fail}", file=sys.stderr)
    """

    stub:
    """
    touch ${assembly_id}.functional_annotations.tsv
    """
}