/*
========================================================================================
    FILTER GENE MODELS
========================================================================================
    Post-BRAKER quality filter to remove clearly artifactual gene predictions.

    Without RNA-seq evidence, BRAKER protein-only mode tends to over-predict
    genes — particularly short single-exon ORFs in repeat-derived regions and
    gene models with impossibly short introns. This module removes predictions
    that fail basic structural plausibility checks while preserving all
    defensible gene models.

    Filters applied (gene is removed if ANY criterion is met):
      1. Single-exon gene with CDS < min_cds_length bp (default: 300 bp / 100 aa)
         → Likely TE ORFs or random predictions
      2. Gene total span < min_gene_length bp (default: 200 bp)
         → Too short to be a real gene
      3. ALL introns in the gene are < min_intron_length bp (default: 10 bp)
         → Structural artifacts from GeneMark

    Removal is at the gene level: if a gene is flagged, all its children
    (mRNA, exon, CDS, intron, start_codon, stop_codon) are also removed.

    Also filters the BRAKER protein FASTA to match, so downstream tools
    (BUSCO, OMArk, functional annotation) operate on the cleaned gene set.

    Inputs:
        - assembly_id:       Sample identifier
        - braker_gff3:       BRAKER3 gene models GFF3
        - braker_proteins:   BRAKER3 predicted proteins FASTA

    Outputs:
        - filtered_gff3:     ${assembly_id}.filtered_genes.gff3
        - filtered_proteins: ${assembly_id}.filtered_proteins.fa
        - removed_gff3:      ${assembly_id}.removed_genes.gff3
        - filter_report:     ${assembly_id}.gene_filter_report.tsv

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process FILTER_GENE_MODELS {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/structural_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(braker_gff3), path(braker_proteins), path(aed_scores)

    output:
    tuple val(assembly_id), path("${assembly_id}.filtered_genes.gff3"),
        emit: filtered_gff3
    tuple val(assembly_id), path("${assembly_id}.filtered_proteins.fa"),
        emit: filtered_proteins
    tuple val(assembly_id), path("${assembly_id}.removed_genes.gff3"),
        emit: removed_gff3
    tuple val(assembly_id), path("${assembly_id}.gene_filter_report.tsv"),
        emit: filter_report
    tuple val(assembly_id), path("${assembly_id}.filtered_aed.tsv"),
        emit: filtered_aed

    script:
    def min_cds_length    = params.filter_min_cds_length   ?: 300
    def min_gene_length   = params.filter_min_gene_length  ?: 200
    def min_intron_length = params.filter_min_intron_length ?: 10
    """
    #!/usr/bin/env python3
    import sys
    from collections import defaultdict, OrderedDict

    MIN_CDS_LENGTH    = ${min_cds_length}
    MIN_GENE_LENGTH   = ${min_gene_length}
    MIN_INTRON_LENGTH = ${min_intron_length}
    ASM = "${assembly_id}"

    # ── Parse GFF3 into gene trees ────────────────────────────────────────────

    def parse_attrs(attr_str):
        d = OrderedDict()
        if not attr_str or attr_str == ".":
            return d
        for token in attr_str.split(";"):
            token = token.strip()
            if "=" in token:
                k, v = token.split("=", 1)
                d[k.strip()] = v.strip()
        return d

    # Store all features grouped by gene
    genes = OrderedDict()
    gene_children = defaultdict(list)
    feature_to_gene = {}
    comments = []

    current_gene_id = None

    with open("${braker_gff3}") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                comments.append(line)
                continue

            cols = line.rstrip("\\n").split("\\t")
            if len(cols) != 9:
                comments.append(line)
                continue

            ftype = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            attrs = parse_attrs(cols[8])
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")

            if ftype == "gene":
                gene_id = fid
                genes[gene_id] = {
                    "line": line,
                    "seqid": cols[0],
                    "start": start,
                    "end": end,
                    "span": end - start + 1,
                    "attrs": attrs,
                }
                feature_to_gene[gene_id] = gene_id
                current_gene_id = gene_id

            elif ftype in ("mRNA", "transcript"):
                gene_id = parent if parent in genes else feature_to_gene.get(parent, current_gene_id)
                feature_to_gene[fid] = gene_id
                gene_children[gene_id].append({
                    "line": line, "type": ftype,
                    "start": start, "end": end, "parent": parent, "id": fid,
                })

            else:
                gene_id = feature_to_gene.get(parent, current_gene_id)
                if fid:
                    feature_to_gene[fid] = gene_id
                gene_children[gene_id].append({
                    "line": line, "type": ftype,
                    "start": start, "end": end, "parent": parent, "id": fid,
                })

    # ── Read AED scores ───────────────────────────────────────────────────────
    import csv
    gene_aed = {}   # gene_id → min AED across transcripts
    tx_aed = {}     # transcript_id → aed

    with open("${aed_scores}") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            gid = row["gene_id"]
            tid = row["transcript_id"]
            aed = float(row["aed"])
            tx_aed[tid] = aed
            if gid not in gene_aed or aed < gene_aed[gid]:
                gene_aed[gid] = aed

    MAX_AED = ${params.filter_max_aed != null ? params.filter_max_aed : 'None'}

    # ── Apply filters ─────────────────────────────────────────────────────────
    removed_genes = {}

    for gene_id, gene_info in genes.items():
        children = gene_children.get(gene_id, [])

        exons = [c for c in children if c["type"] == "exon"]
        cds_features = [c for c in children if c["type"] == "CDS"]
        introns = [c for c in children if c["type"] == "intron"]
        mrnas = [c for c in children if c["type"] in ("mRNA", "transcript")]

        gene_span = gene_info["span"]

        # Count exons per mRNA
        exons_per_mrna = defaultdict(int)
        for e in exons:
            exons_per_mrna[e["parent"]] += 1

        total_cds_bp = sum(c["end"] - c["start"] + 1 for c in cds_features)

        all_single_exon = all(
            exons_per_mrna.get(m["id"], 1) == 1 for m in mrnas
        ) if mrnas else (len(exons) <= 1)

        intron_lengths = [c["end"] - c["start"] + 1 for c in introns]

        # Filter 1: Short single-exon gene
        if all_single_exon and total_cds_bp < MIN_CDS_LENGTH and total_cds_bp > 0:
            removed_genes[gene_id] = f"single_exon_short_cds({total_cds_bp}bp<{MIN_CDS_LENGTH}bp)"
            continue

        # Filter 2: Very short gene span
        if gene_span < MIN_GENE_LENGTH:
            removed_genes[gene_id] = f"short_gene_span({gene_span}bp<{MIN_GENE_LENGTH}bp)"
            continue

        # Filter 3: All introns impossibly short
        if intron_lengths and all(il < MIN_INTRON_LENGTH for il in intron_lengths):
            max_il = max(intron_lengths)
            removed_genes[gene_id] = f"all_introns_short(max={max_il}bp<{MIN_INTRON_LENGTH}bp)"
            continue
        
        # Filter 4: AED threshold (if configured)
        if MAX_AED is not None:
            best = gene_aed.get(gene_id, 1.0)
            if best > MAX_AED:
                removed_genes[gene_id] = f"high_AED({best:.3f}>{MAX_AED})"
                continue

    # ── Collect kept mRNA/transcript IDs for protein filtering ────────────────
    kept_mrna_ids = set()
    removed_mrna_ids = set()

    for gene_id in genes:
        children = gene_children.get(gene_id, [])
        mrna_ids = [c["id"] for c in children if c["type"] in ("mRNA", "transcript")]

        if gene_id in removed_genes:
            removed_mrna_ids.update(mrna_ids)
        else:
            kept_mrna_ids.update(mrna_ids)

    # Also keep gene IDs themselves (BRAKER protein headers can use either)
    kept_gene_ids = set(genes.keys()) - set(removed_genes.keys())
    kept_ids = kept_mrna_ids | kept_gene_ids

    # ── Write GFF3 outputs ────────────────────────────────────────────────────
    with open(f"{ASM}.filtered_genes.gff3", "w") as kept_f, \\
         open(f"{ASM}.removed_genes.gff3", "w") as removed_f:

        for c in comments:
            kept_f.write(c)

        for gene_id, gene_info in genes.items():
            target = removed_f if gene_id in removed_genes else kept_f
            target.write(gene_info["line"])
            for child in gene_children.get(gene_id, []):
                target.write(child["line"])

    # ── Filter protein FASTA ──────────────────────────────────────────────────
    #   BRAKER protein headers typically use mRNA ID (e.g. g1.t1) or gene ID.
    #   We keep a protein if its header ID matches any kept mRNA or gene ID.
    n_prot_kept = 0
    n_prot_removed = 0

    with open(f"{ASM}.filtered_proteins.fa", "w") as out:
        keep_current = False
        with open("${braker_proteins}") as fh:
            for line in fh:
                if line.startswith(">"):
                    prot_id = line.strip().split()[0][1:]
                    # Check against kept IDs (try exact match first, then base)
                    if prot_id in kept_ids:
                        keep_current = True
                    else:
                        # Try stripping .t1 suffix for gene-level match
                        base = prot_id.rsplit(".", 1)[0] if "." in prot_id else prot_id
                        keep_current = base in kept_ids or prot_id in kept_ids

                    if keep_current:
                        n_prot_kept += 1
                    else:
                        n_prot_removed += 1

                if keep_current:
                    out.write(line)

    # ── Write filter report ───────────────────────────────────────────────────
    with open(f"{ASM}.gene_filter_report.tsv", "w") as out:
        out.write("gene_id\\tseqid\\tstart\\tend\\tspan\\tn_exons\\tcds_bp\\t"
                  "n_introns\\tmin_intron\\tstatus\\treason\\n")

        for gene_id, gene_info in genes.items():
            children = gene_children.get(gene_id, [])
            exons = [c for c in children if c["type"] == "exon"]
            cds_features = [c for c in children if c["type"] == "CDS"]
            introns = [c for c in children if c["type"] == "intron"]

            total_cds = sum(c["end"] - c["start"] + 1 for c in cds_features)
            intron_lens = [c["end"] - c["start"] + 1 for c in introns]
            min_intron = min(intron_lens) if intron_lens else 0

            status = "removed" if gene_id in removed_genes else "kept"
            reason = removed_genes.get(gene_id, "-")

            out.write(f"{gene_id}\\t{gene_info['seqid']}\\t{gene_info['start']}\\t"
                      f"{gene_info['end']}\\t{gene_info['span']}\\t{len(exons)}\\t"
                      f"{total_cds}\\t{len(introns)}\\t{min_intron}\\t{status}\\t{reason}\\n")

    # ── Write filtered AED scores ────────────────────────────────────────────
    with open(f"{ASM}.filtered_aed.tsv", "w") as out:
        out.write("gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\t"
                "protein_identity\\tn_evidence_hits\\n")
        with open("${aed_scores}") as fh:
            reader = csv.DictReader(fh, delimiter="\\t")
            for row in reader:
                if row["gene_id"] not in removed_genes:
                    out.write("\\t".join(row[c] for c in reader.fieldnames) + "\\n")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_total = len(genes)
    n_removed = len(removed_genes)
    n_kept = n_total - n_removed

    reasons = defaultdict(int)
    for r in removed_genes.values():
        key = r.split("(")[0]
        reasons[key] += 1

    print(f"Gene model filter: {n_total} total → {n_kept} kept, {n_removed} removed",
          file=sys.stderr)
    for reason, count in sorted(reasons.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count}", file=sys.stderr)
    print(f"Proteins: {n_prot_kept} kept, {n_prot_removed} removed", file=sys.stderr)
    """

    stub:
    """
    touch ${assembly_id}.filtered_genes.gff3
    touch ${assembly_id}.filtered_proteins.fa
    touch ${assembly_id}.removed_genes.gff3
    touch ${assembly_id}.gene_filter_report.tsv
    printf "gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\tprotein_identity\\tn_evidence_hits\\n" \\
        > ${assembly_id}.filtered_aed.tsv
    """
}