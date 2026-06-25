process RENAME_GENE_MODELS {
    tag "${assembly_id}"
    label 'rename_gene_models'
    publishDir "${params.outdir}/intermediate/renamed_gene_models", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id),
          path(flagged_gff3),
          path(proteins_fasta),
          path(aed_scores)

    output:
    tuple val(assembly_id),
          path("${assembly_id}.renamed.gff3"),
          emit: renamed_gff3

    tuple val(assembly_id),
          path("${assembly_id}.renamed.proteins.faa"),
          emit: renamed_proteins

    tuple val(assembly_id),
          path("${assembly_id}.renamed.aed.tsv"),
          emit: renamed_aed

    tuple val(assembly_id),
          path("${assembly_id}.gene_model_id_map.tsv"),
          emit: id_map

    script:
    """\
    #!/usr/bin/env python3

    import csv
    import re
    import sys
    import textwrap
    from collections import defaultdict, OrderedDict, Counter

    ASM = "${assembly_id}"

    gff_in = "${flagged_gff3}"
    prot_in = "${proteins_fasta}"
    aed_in = "${aed_scores}"

    gff_out = f"{ASM}.renamed.gff3"
    prot_out = f"{ASM}.renamed.proteins.faa"
    aed_out = f"{ASM}.renamed.aed.tsv"
    map_out = f"{ASM}.gene_model_id_map.tsv"

    def parse_attributes(attr_str):
        attrs = OrderedDict()
        if not attr_str or attr_str == ".":
            return attrs

        for token in attr_str.split(";"):
            token = token.strip()
            if not token:
                continue

            if "=" in token:
                key, value = token.split("=", 1)
                attrs[key] = value
            else:
                attrs[token] = ""

        return attrs

    def format_attributes(attrs):
        if not attrs:
            return "."

        parts = []
        for key, value in attrs.items():
            if value == "":
                parts.append(key)
            else:
                parts.append(f"{key}={value}")

        return ";".join(parts)

    def clean_feature_type(ftype):
        return re.sub(r"[^A-Za-z0-9_]", "_", ftype)

    # -------------------------------------------------------------------------
    # Read GFF
    # -------------------------------------------------------------------------

    raw_records = []
    gff_headers = []

    with open(gff_in) as fh:
        for line in fh:
            line = line.rstrip("\\n")

            if not line:
                continue

            if line.startswith("#"):
                gff_headers.append(line)
                continue

            cols = line.split("\\t")

            if len(cols) != 9:
                continue

            attrs = parse_attributes(cols[8])
            raw_records.append([cols[:8], attrs])

    # -------------------------------------------------------------------------
    # Rename IDs
    # -------------------------------------------------------------------------

    id_map = {}
    gene_map = {}
    mrna_map = {}

    map_rows = []

    gene_counter = 0
    mrna_counters = defaultdict(int)
    child_counters = defaultdict(int)

    # Pass 1: rename gene features
    for cols8, attrs in raw_records:
        ftype = cols8[2]
        old_id = attrs.get("ID", "")

        if ftype != "gene":
            continue

        if not old_id:
            raise SystemExit(
                f"ERROR: gene feature without ID at {cols8[0]}:{cols8[3]}-{cols8[4]}"
            )

        gene_counter += 1
        new_id = f"{ASM}_g{gene_counter}"

        id_map[old_id] = new_id
        gene_map[old_id] = new_id

        attrs["ID"] = new_id

        map_rows.append({
            "assembly_id": ASM,
            "feature_type": ftype,
            "old_id": old_id,
            "new_id": new_id,
            "old_parent": ".",
            "new_parent": ".",
            "seqid": cols8[0],
            "start": cols8[3],
            "end": cols8[4],
            "strand": cols8[6],
        })

    # Pass 2: rename mRNA/transcript features
    for cols8, attrs in raw_records:
        ftype = cols8[2]
        old_id = attrs.get("ID", "")
        old_parent = attrs.get("Parent", "")

        if ftype not in ("mRNA", "transcript"):
            continue

        if not old_id:
            raise SystemExit(
                f"ERROR: {ftype} feature without ID at {cols8[0]}:{cols8[3]}-{cols8[4]}"
            )

        if not old_parent:
            raise SystemExit(f"ERROR: {ftype} {old_id} has no Parent")

        first_old_parent = old_parent.split(",")[0]

        if first_old_parent not in gene_map:
            raise SystemExit(
                f"ERROR: {ftype} {old_id} has Parent={first_old_parent}, "
                "but that gene was not renamed"
            )

        new_parent = gene_map[first_old_parent]

        mrna_counters[new_parent] += 1
        new_id = f"{new_parent}.t{mrna_counters[new_parent]}"

        id_map[old_id] = new_id
        mrna_map[old_id] = new_id

        attrs["ID"] = new_id
        attrs["Parent"] = new_parent

        map_rows.append({
            "assembly_id": ASM,
            "feature_type": ftype,
            "old_id": old_id,
            "new_id": new_id,
            "old_parent": first_old_parent,
            "new_parent": new_parent,
            "seqid": cols8[0],
            "start": cols8[3],
            "end": cols8[4],
            "strand": cols8[6],
        })

    # Pass 3: update child feature Parents and rename child IDs when present
    for cols8, attrs in raw_records:
        ftype = cols8[2]

        if ftype in ("gene", "mRNA", "transcript"):
            continue

        old_parent_raw = attrs.get("Parent", "")
        old_id = attrs.get("ID", "")

        new_parent_raw = old_parent_raw

        if old_parent_raw:
            old_parents = old_parent_raw.split(",")
            new_parents = [id_map.get(parent, parent) for parent in old_parents]
            new_parent_raw = ",".join(new_parents)
            attrs["Parent"] = new_parent_raw

        if old_id:
            first_new_parent = new_parent_raw.split(",")[0] if new_parent_raw else "unparented"
            safe_ftype = clean_feature_type(ftype)

            child_counters[(first_new_parent, safe_ftype)] += 1
            new_id = f"{first_new_parent}.{safe_ftype}{child_counters[(first_new_parent, safe_ftype)]}"

            id_map[old_id] = new_id
            attrs["ID"] = new_id

            map_rows.append({
                "assembly_id": ASM,
                "feature_type": ftype,
                "old_id": old_id,
                "new_id": new_id,
                "old_parent": old_parent_raw if old_parent_raw else ".",
                "new_parent": new_parent_raw if new_parent_raw else ".",
                "seqid": cols8[0],
                "start": cols8[3],
                "end": cols8[4],
                "strand": cols8[6],
            })

    # -------------------------------------------------------------------------
    # Validate renamed IDs
    # -------------------------------------------------------------------------

    new_ids = []

    for cols8, attrs in raw_records:
        feature_id = attrs.get("ID", "")
        if feature_id:
            new_ids.append(feature_id)

    duplicated = [feature_id for feature_id, n in Counter(new_ids).items() if n > 1]

    if duplicated:
        raise SystemExit(
            "ERROR: duplicate renamed IDs detected. First duplicates: "
            + ", ".join(duplicated[:20])
        )

    n_genes = sum(1 for cols8, attrs in raw_records if cols8[2] == "gene")
    n_mrnas = sum(1 for cols8, attrs in raw_records if cols8[2] in ("mRNA", "transcript"))

    if n_genes != gene_counter:
        raise SystemExit(
            f"ERROR: expected to rename {n_genes} genes, but renamed {gene_counter}"
        )

    if n_mrnas != len(mrna_map):
        raise SystemExit(
            f"ERROR: expected to rename {n_mrnas} mRNAs/transcripts, "
            f"but renamed {len(mrna_map)}"
        )

    # -------------------------------------------------------------------------
    # Write renamed GFF3
    # -------------------------------------------------------------------------

    with open(gff_out, "w") as out:
        wrote_version = False

        for header in gff_headers:
            if header.startswith("##gff-version"):
                out.write(header + "\\n")
                wrote_version = True
                break

        if not wrote_version:
            out.write("##gff-version 3\\n")

        out.write(f"# Protein-coding gene models renamed by RENAME_GENE_MODELS for {ASM}\\n")
        out.write(f"# Protein-coding genes: {gene_counter}\\n")
        out.write(f"# Protein-coding mRNAs/transcripts: {len(mrna_map)}\\n")

        for cols8, attrs in raw_records:
            out.write("\\t".join(cols8 + [format_attributes(attrs)]) + "\\n")

    # -------------------------------------------------------------------------
    # Write ID map
    # -------------------------------------------------------------------------

    map_cols = [
        "assembly_id",
        "feature_type",
        "old_id",
        "new_id",
        "old_parent",
        "new_parent",
        "seqid",
        "start",
        "end",
        "strand",
    ]

    with open(map_out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=map_cols, delimiter="\\t")
        writer.writeheader()

        for row in map_rows:
            writer.writerow(row)

    # -------------------------------------------------------------------------
    # Rename protein FASTA headers
    # -------------------------------------------------------------------------

    n_proteins = 0
    unmatched_proteins = []

    current_old = None
    current_new = None
    seq_parts = []

    def write_fasta_record(out, header, sequence):
        out.write(f">{header}\\n")
        for chunk in textwrap.wrap(sequence, 60):
            out.write(chunk + "\\n")

    with open(prot_in) as fh, open(prot_out, "w") as out:
        for line in fh:
            line = line.rstrip("\\n")

            if line.startswith(">"):
                if current_old is not None:
                    write_fasta_record(out, current_new, "".join(seq_parts))
                    n_proteins += 1

                current_old = line[1:].split()[0]
                current_new = mrna_map.get(current_old)

                if current_new is None:
                    unmatched_proteins.append(current_old)
                    current_new = current_old

                seq_parts = []

            else:
                seq_parts.append(line.strip())

        if current_old is not None:
            write_fasta_record(out, current_new, "".join(seq_parts))
            n_proteins += 1

    if unmatched_proteins:
        raise SystemExit(
            "ERROR: protein FASTA contains IDs not present as mRNA/transcript IDs in the GFF. "
            f"Count={len(unmatched_proteins)}. First examples: "
            + ", ".join(unmatched_proteins[:20])
        )

    if n_proteins != len(mrna_map):
        raise SystemExit(
            f"ERROR: renamed protein FASTA has {n_proteins} records, "
            f"but renamed GFF has {len(mrna_map)} mRNAs/transcripts"
        )

    # -------------------------------------------------------------------------
    # Rename AED transcript IDs
    # -------------------------------------------------------------------------

    n_aed_written = 0
    n_aed_skipped = 0

    with open(aed_in) as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        fieldnames = reader.fieldnames

        if fieldnames is None:
            fieldnames = ["transcript_id"]

        out_fields = list(fieldnames)

        if "original_transcript_id" not in out_fields:
            out_fields.append("original_transcript_id")

        with open(aed_out, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=out_fields, delimiter="\\t")
            writer.writeheader()

            for row in reader:
                old_tid = row.get("transcript_id", "").strip()

                if not old_tid:
                    continue

                new_tid = mrna_map.get(old_tid)

                if new_tid is None:
                    n_aed_skipped += 1
                    continue

                row["original_transcript_id"] = old_tid
                row["transcript_id"] = new_tid

                writer.writerow(row)
                n_aed_written += 1

    print(
        f"Renamed {gene_counter} genes, {len(mrna_map)} mRNAs/proteins, "
        f"{n_aed_written} AED rows for {ASM}. "
        f"Skipped {n_aed_skipped} AED rows not present after filtering.",
        file=sys.stderr
    )
    """.stripIndent()

    stub:
    """
    touch ${assembly_id}.renamed.gff3
    touch ${assembly_id}.renamed.proteins.faa
    touch ${assembly_id}.renamed.aed.tsv
    touch ${assembly_id}.gene_model_id_map.tsv
    """
}