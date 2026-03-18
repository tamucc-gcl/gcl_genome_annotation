/*
========================================================================================
    ANNOTATE GFF3 MODULE
========================================================================================
    Merge structural (BRAKER), tRNAscan, Infernal, telomere, centromere, and
    repeat GFF3 annotations into a single unified GFF3 with:
      - Consistent ID namespacing ({assembly_id}_g{N}, _tRNA{N}, _ncRNA{N})
      - Functional attributes (Dbxref, product, Ontology_term, Note) on
        protein-coding mRNA features from the merged functional table
      - Product name hierarchy:
            1. DIAMOND SwissProt hit description
            2. eggNOG preferred_name + description
            3. KofamScan KO definition
            4. "hypothetical protein"

    Inputs:
        - assembly_id:          Sample identifier
        - braker_gff3:          BRAKER3 gene models GFF3
        - trnascan_gff3:        tRNAscan-SE predictions GFF3
        - infernal_gff3:        Infernal/Rfam ncRNA predictions GFF3
        - telomere_gff3:        tidk telomere annotations GFF3
        - centromere_gff3:      centroAnno centromere annotations GFF3
        - repeat_gff:           RepeatMasker repeat annotations GFF
        - functional_tsv:       Merged functional annotations TSV
                                (from MERGE_FUNCTIONAL_ANNOTATIONS)

    Outputs:
        - annotated_gff3:  ${assembly_id}.annotated.gff3
        - gene_summary:    ${assembly_id}.gene_product_summary.tsv

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process ANNOTATE_GFF3 {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id),
          path(braker_gff3),
          path(trnascan_gff3),
          path(infernal_gff3),
          path(telomere_gff3),
          path(centromere_gff3),
          path(repeat_gff),
          path(functional_tsv)

    output:
    tuple val(assembly_id),
          path("${assembly_id}.annotated.gff3"),
          emit: annotated_gff3
    tuple val(assembly_id),
          path("${assembly_id}.gene_product_summary.tsv"),
          emit: gene_summary

    script:
    """
    #!/usr/bin/env python3
    \"\"\"
    Merge BRAKER, tRNAscan, Infernal, telomere, centromere, and repeat GFF3s
    into a single annotated GFF3.

    1. Parse & renumber BRAKER features  → {asm}_g{N}, mRNA → {asm}_g{N}.t{M}
    2. Parse & renumber tRNAscan features → {asm}_tRNA{N}
    3. Parse & renumber Infernal features → {asm}_ncRNA{N}
    4. Parse telomere GFF3               → {asm}_telo{N}
    5. Parse centromere GFF3             → pass through (already has IDs)
    6. Parse repeat GFF                  → pass through (RepeatMasker features)
    7. Load functional table, decorate mRNA features with product/Dbxref/GO/Note
    8. Propagate best product name to parent gene feature
    9. Write merged GFF3 sorted by seqid, start
    \"\"\"
    import csv, re, sys, urllib.parse
    from collections import defaultdict, OrderedDict

    ASM = "${assembly_id}"

    # ═══════════════════════════════════════════════════════════════════════════
    # Helpers
    # ═══════════════════════════════════════════════════════════════════════════

    def parse_attributes(attr_str):
        \"\"\"Parse GFF3 column-9 into an OrderedDict (preserves key order).\"\"\"
        d = OrderedDict()
        if not attr_str or attr_str == ".":
            return d
        for token in attr_str.split(";"):
            token = token.strip()
            if not token:
                continue
            if "=" in token:
                k, v = token.split("=", 1)
                d[k] = v
            else:
                d[token] = ""
        return d

    def format_attributes(d):
        \"\"\"Serialize an OrderedDict back to GFF3 column-9 string.\"\"\"
        parts = []
        for k, v in d.items():
            if v == "":
                parts.append(k)
            else:
                parts.append(f"{k}={v}")
        return ";".join(parts) if parts else "."

    def gff_sort_key(rec):
        \"\"\"Sort GFF records by seqid then start position.\"\"\"
        return (rec[0], int(rec[3]))

    # ═══════════════════════════════════════════════════════════════════════════
    # 1. Parse BRAKER GFF3 — renumber genes and propagate new IDs to children
    # ═══════════════════════════════════════════════════════════════════════════

    braker_records = []
    braker_id_map = {}
    braker_gene_counter = 0
    braker_mrna_counters = defaultdict(int)

    raw_braker = []
    with open("${braker_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) != 9:
                continue
            attrs = parse_attributes(cols[8])
            raw_braker.append((cols[:8], attrs))

    # Pass 1a: assign new gene IDs
    for cols8, attrs in raw_braker:
        ftype = cols8[2]
        old_id = attrs.get("ID", "")
        if ftype == "gene":
            braker_gene_counter += 1
            new_id = f"{ASM}_g{braker_gene_counter}"
            braker_id_map[old_id] = new_id
            attrs["ID"] = new_id

    # Pass 1b: assign new mRNA/transcript IDs and update Parent refs
    for cols8, attrs in raw_braker:
        ftype = cols8[2]
        old_id = attrs.get("ID", "")
        parent = attrs.get("Parent", "")
        if ftype in ("mRNA", "transcript"):
            new_parent = braker_id_map.get(parent, parent)
            braker_mrna_counters[new_parent] += 1
            new_id = f"{new_parent}.t{braker_mrna_counters[new_parent]}"
            braker_id_map[old_id] = new_id
            attrs["ID"] = new_id
            attrs["Parent"] = new_parent

    # Pass 1c: update all child features (exon, CDS, etc.) Parent references
    for cols8, attrs in raw_braker:
        ftype = cols8[2]
        if ftype in ("gene", "mRNA", "transcript"):
            continue
        parent = attrs.get("Parent", "")
        if parent in braker_id_map:
            attrs["Parent"] = braker_id_map[parent]
        old_id = attrs.get("ID", "")
        if old_id:
            new_parent = attrs.get("Parent", parent)
            attrs["ID"] = f"{new_parent}.{ftype}{old_id[-1] if old_id[-1].isdigit() else '1'}"

    for cols8, attrs in raw_braker:
        braker_records.append(cols8 + [attrs])

    # Build mRNA→gene parent mapping for product propagation later
    mrna_to_gene = {}
    for cols8, attrs in raw_braker:
        ftype = cols8[2]
        if ftype in ("mRNA", "transcript"):
            mrna_to_gene[attrs["ID"]] = attrs.get("Parent", "")

    old_mrna_to_new = {}
    for old_id, new_id in braker_id_map.items():
        if any(new_id.endswith(f".t{i}") for i in range(1, 100)):
            old_mrna_to_new[old_id] = new_id

    # ═══════════════════════════════════════════════════════════════════════════
    # 2. Parse tRNAscan GFF3 — renumber
    # ═══════════════════════════════════════════════════════════════════════════

    trna_records = []
    trna_counter = 0
    trna_pseudo_counter = 0
    trna_id_map = {}

    raw_trna = []
    with open("${trnascan_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) != 9:
                continue
            attrs = parse_attributes(cols[8])
            raw_trna.append((cols[:8], attrs))

    # Pass 1: assign new IDs to top-level features (tRNA, pseudogene, or gene)
    for cols8, attrs in raw_trna:
        ftype = cols8[2]
        old_id = attrs.get("ID", "")
        parent = attrs.get("Parent", "")

        if ftype == "gene":
            # Some tRNAscan versions emit gene wrappers
            trna_counter += 1
            new_id = f"{ASM}_tRNA{trna_counter}"
            trna_id_map[old_id] = new_id
            attrs["ID"] = new_id

        elif ftype == "tRNA":
            if parent and parent in trna_id_map:
                # tRNA is child of a gene wrapper
                new_parent = trna_id_map[parent]
                new_id = f"{new_parent}.trna1"
                trna_id_map[old_id] = new_id
                attrs["ID"] = new_id
                attrs["Parent"] = new_parent
            else:
                # tRNA is top-level (no gene parent) — common in tRNAscan v2.x
                trna_counter += 1
                new_id = f"{ASM}_tRNA{trna_counter}"
                trna_id_map[old_id] = new_id
                attrs["ID"] = new_id
                # Remove empty/dangling Parent if present
                if "Parent" in attrs and not attrs["Parent"]:
                    del attrs["Parent"]

        elif ftype == "pseudogene":
            trna_pseudo_counter += 1
            new_id = f"{ASM}_tRNA_pseudo{trna_pseudo_counter}"
            trna_id_map[old_id] = new_id
            attrs["ID"] = new_id

    # Pass 2: update child features (exons) to point to new parent IDs
    for cols8, attrs in raw_trna:
        ftype = cols8[2]
        if ftype in ("gene", "tRNA", "pseudogene"):
            continue
        old_parent = attrs.get("Parent", "")
        if old_parent in trna_id_map:
            new_parent = trna_id_map[old_parent]
            attrs["Parent"] = new_parent
            old_id = attrs.get("ID", "")
            if old_id:
                exon_match = re.search(r'exon(\\d+)', old_id)
                exon_num = exon_match.group(1) if exon_match else "1"
                attrs["ID"] = f"{new_parent}.exon{exon_num}"

    for cols8, attrs in raw_trna:
        trna_records.append(cols8 + [attrs])

    # ═══════════════════════════════════════════════════════════════════════════
    # 3. Parse Infernal GFF3 — renumber
    # ═══════════════════════════════════════════════════════════════════════════

    ncrna_records = []
    ncrna_counter = 0

    with open("${infernal_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            attrs = parse_attributes(cols[8])
            ncrna_counter += 1

            new_id = f"{ASM}_ncRNA{ncrna_counter}"
            attrs["ID"] = new_id

            rfam_model = cols[2]
            attrs["Name"] = rfam_model
            cols[2] = "ncRNA"

            ncrna_records.append(cols[:8] + [attrs])

    # ═══════════════════════════════════════════════════════════════════════════
    # 3b. Parse telomere GFF3
    # ═══════════════════════════════════════════════════════════════════════════

    telo_records = []
    telo_counter = 0

    with open("${telomere_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            attrs = parse_attributes(cols[8])
            telo_counter += 1

            new_id = f"{ASM}_telo{telo_counter}"
            attrs["ID"] = new_id

            telo_records.append(cols[:8] + [attrs])

    # ═══════════════════════════════════════════════════════════════════════════
    # 3c. Parse centromere GFF3
    # ═══════════════════════════════════════════════════════════════════════════

    centro_records = []
    centro_counter = 0

    with open("${centromere_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            attrs = parse_attributes(cols[8])
            centro_counter += 1
            # Keep existing IDs from CENTROMERE_TO_GFF3
            centro_records.append(cols[:8] + [attrs])

    # ═══════════════════════════════════════════════════════════════════════════
    # 3d. Parse repeat GFF (RepeatMasker output)
    # ═══════════════════════════════════════════════════════════════════════════

    repeat_records = []
    repeat_counter = 0

    with open("${repeat_gff}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            attrs = parse_attributes(cols[8])
            repeat_counter += 1

            # Assign assembly-namespaced IDs
            new_id = f"{ASM}_repeat{repeat_counter}"
            attrs["ID"] = new_id

            repeat_records.append(cols[:8] + [attrs])

    # ═══════════════════════════════════════════════════════════════════════════
    # 4. Load functional annotations
    # ═══════════════════════════════════════════════════════════════════════════

    func = {}
    with open("${functional_tsv}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            pid = row.get("protein_id", "").strip()
            if pid:
                func[pid] = row

    def get_product_name(rec):
        \"\"\"Determine product name using hierarchy:
           1. DIAMOND SwissProt description (clean up OS=... suffix)
           2. eggNOG preferred_name + description
           3. KofamScan KO definition
           4. hypothetical protein
        \"\"\"
        d_desc = rec.get("diamond_desc", "-").strip()
        if d_desc and d_desc != "-":
            cleaned = re.sub(r"\\s+OS=.*", "", d_desc)
            cleaned = re.sub(r"^\\S+\\s+", "", cleaned, count=1).strip()
            if cleaned:
                return cleaned

        e_name = rec.get("eggnog_name", "-").strip()
        e_desc = rec.get("eggnog_desc", "-").strip()
        if e_name and e_name != "-" and e_desc and e_desc != "-":
            return f"{e_name}: {e_desc}"
        if e_desc and e_desc != "-":
            return e_desc
        if e_name and e_name != "-":
            return e_name

        k_def = rec.get("kofam_def", "-").strip()
        if k_def and k_def != "-":
            return k_def

        return "hypothetical protein"

    def collect_dbxrefs(rec):
        \"\"\"Build Dbxref list from available annotations.\"\"\"
        xrefs = []
        d_hit = rec.get("diamond_hit", "-").strip()
        if d_hit and d_hit != "-":
            xrefs.append(f"UniProtKB/Swiss-Prot:{d_hit}")

        e_ogs = rec.get("eggnog_OGs", "-").strip()
        if e_ogs and e_ogs != "-":
            for og in e_ogs.split(","):
                og = og.strip().split("@")[0]
                if og:
                    xrefs.append(f"eggNOG:{og}")

        e_ko = rec.get("eggnog_KEGG_ko", "-").strip()
        if e_ko and e_ko != "-":
            for ko in e_ko.split(","):
                ko = ko.strip()
                if ko.startswith("ko:"):
                    xrefs.append(f"KEGG:{ko}")

        k_ko = rec.get("kofam_KO", "-").strip()
        if k_ko and k_ko != "-":
            xrefs.append(f"KEGG_KO:{k_ko}")

        ipr_acc = rec.get("ipr_accessions", "-").strip()
        if ipr_acc and ipr_acc != "-":
            for acc in ipr_acc.split(";"):
                acc = acc.strip()
                if acc and acc != "-":
                    xrefs.append(f"InterPro:{acc}")

        return xrefs

    def collect_go_terms(rec):
        \"\"\"Collect GO terms from InterProScan and eggNOG (deduplicated).\"\"\"
        terms = set()
        for field in ("ipr_go", "eggnog_GO"):
            raw = rec.get(field, "-").strip()
            if raw and raw != "-":
                for t in re.split(r"[;,|]+", raw):
                    t = t.strip()
                    if t.startswith("GO:"):
                        terms.add(t)
        return sorted(terms)

    def build_note(rec):
        \"\"\"Build a Note attribute with supplementary info.\"\"\"
        parts = []
        ds_type = rec.get("deepsig_type", "-").strip()
        if ds_type and ds_type != "-" and ds_type.lower() != "other":
            ds_start = rec.get("deepsig_start", "?")
            ds_end = rec.get("deepsig_end", "?")
            parts.append(f"Signal peptide ({ds_type}) {ds_start}-{ds_end}")

        cazy = rec.get("eggnog_CAZy", "-").strip()
        if cazy and cazy != "-":
            parts.append(f"CAZy:{cazy}")

        kpath = rec.get("eggnog_KEGG_path", "-").strip()
        if kpath and kpath != "-":
            parts.append(f"KEGG_pathway:{kpath}")

        return "; ".join(parts) if parts else None

    # ═══════════════════════════════════════════════════════════════════════════
    # 5. Decorate mRNA features with functional attributes
    # ═══════════════════════════════════════════════════════════════════════════

    gene_products = {}

    for rec_list in [braker_records]:
        for entry in rec_list:
            cols8 = entry[:8]
            attrs = entry[8]
            ftype = cols8[2]

            if ftype not in ("mRNA", "transcript"):
                continue

            new_mrna_id = attrs["ID"]

            func_rec = None
            for old_id, new_id in braker_id_map.items():
                if new_id == new_mrna_id:
                    if old_id in func:
                        func_rec = func[old_id]
                        break
                    for suffix in ("", ".t1", "_t1"):
                        candidate = old_id + suffix
                        if candidate in func:
                            func_rec = func[candidate]
                            break
                    base = old_id.rsplit(".", 1)[0] if "." in old_id else old_id
                    if base in func:
                        func_rec = func[base]
                    break

            if not func_rec:
                continue

            product = get_product_name(func_rec)
            attrs["product"] = urllib.parse.quote(product, safe=" :/()-,.'")

            xrefs = collect_dbxrefs(func_rec)
            if xrefs:
                existing = attrs.get("Dbxref", "")
                all_xrefs = ([existing] if existing else []) + xrefs
                attrs["Dbxref"] = ",".join(all_xrefs)

            go_terms = collect_go_terms(func_rec)
            if go_terms:
                attrs["Ontology_term"] = ",".join(go_terms)

            note = build_note(func_rec)
            if note:
                attrs["Note"] = urllib.parse.quote(note, safe=" :/()-,.'")

            parent_gene = mrna_to_gene.get(new_mrna_id, "")
            if parent_gene:
                current = gene_products.get(parent_gene, "hypothetical protein")
                if current == "hypothetical protein" or product != "hypothetical protein":
                    gene_products[parent_gene] = product

    # Propagate product to gene features
    for entry in braker_records:
        cols8 = entry[:8]
        attrs = entry[8]
        if cols8[2] == "gene":
            gene_id = attrs["ID"]
            if gene_id in gene_products:
                attrs["product"] = urllib.parse.quote(
                    gene_products[gene_id], safe=" :/()-,.'")

    # ═══════════════════════════════════════════════════════════════════════════
    # 6. Write merged GFF3
    # ═══════════════════════════════════════════════════════════════════════════

    all_records = []

    for entry in braker_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        cols8[1] = "BRAKER3"
        all_records.append(cols8 + [format_attributes(attrs)])

    for entry in trna_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        cols8[1] = "tRNAscan-SE"
        all_records.append(cols8 + [format_attributes(attrs)])

    for entry in ncrna_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        cols8[1] = "Infernal"
        all_records.append(cols8 + [format_attributes(attrs)])

    for entry in telo_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        all_records.append(cols8 + [format_attributes(attrs)])

    for entry in centro_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        all_records.append(cols8 + [format_attributes(attrs)])

    for entry in repeat_records:
        cols8 = list(entry[:8])
        attrs = entry[8]
        cols8[1] = "RepeatMasker"
        all_records.append(cols8 + [format_attributes(attrs)])

    # Sort by seqid, then start position
    all_records.sort(key=lambda r: (r[0], int(r[3])))

    with open(f"{ASM}.annotated.gff3", "w") as out:
        out.write("##gff-version 3\\n")
        out.write(f"# Genome annotation for {ASM}\\n")
        out.write(f"# Sources: BRAKER3 (protein-coding), tRNAscan-SE (tRNA), "
                  f"Infernal/Rfam (ncRNA), tidk (telomeres), "
                  f"centroAnno (centromeres), RepeatMasker (repeats)\\n")
        out.write(f"# Protein-coding genes: {braker_gene_counter}\\n")
        out.write(f"# tRNA genes: {trna_counter}\\n")
        out.write(f"# ncRNA features: {ncrna_counter}\\n")
        out.write(f"# Telomeric regions: {telo_counter}\\n")
        out.write(f"# Centromere features: {centro_counter}\\n")
        out.write(f"# Repeat features: {repeat_counter}\\n")

        current_seqid = None
        for rec in all_records:
            if rec[0] != current_seqid:
                if current_seqid is not None:
                    out.write("###\\n")
                current_seqid = rec[0]
            out.write("\\t".join(str(c) for c in rec) + "\\n")
        out.write("###\\n")

    # ═══════════════════════════════════════════════════════════════════════════
    # 7. Write gene product summary TSV
    # ═══════════════════════════════════════════════════════════════════════════

    with open(f"{ASM}.gene_product_summary.tsv", "w") as out:
        out.write("gene_id\\tproduct\\n")
        for entry in braker_records:
            attrs = entry[8]
            if entry[2] == "gene":
                gene_id = attrs["ID"]
                product = urllib.parse.unquote(attrs.get("product", "hypothetical protein"))
                out.write(f"{gene_id}\\t{product}\\n")

    # ── Summary stats ────────────────────────────────────────────────────────
    n_annotated = sum(1 for p in gene_products.values() if p != "hypothetical protein")
    n_hypo = braker_gene_counter - n_annotated
    print(f"Merged GFF3: {braker_gene_counter} protein-coding genes, "
          f"{trna_counter} tRNAs, {ncrna_counter} ncRNAs, "
          f"{telo_counter} telomeric regions, "
          f"{centro_counter} centromere features, "
          f"{repeat_counter} repeat features", file=sys.stderr)
    print(f"Functional annotation: {n_annotated} named, {n_hypo} hypothetical", file=sys.stderr)
    """

    stub:
    """
    touch ${assembly_id}.annotated.gff3
    touch ${assembly_id}.gene_product_summary.tsv
    """
}