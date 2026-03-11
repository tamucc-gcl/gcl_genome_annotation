/*
========================================================================================
    centroAnno MODULE
========================================================================================
    De novo centromere and tandem repeat annotation using centroAnno

    centroAnno identifies centromeric satellite arrays and their higher-order
    repeat (HOR) structure directly from genome assemblies without requiring
    prior knowledge of monomer templates.

    Mode: anno-asm (whole-chromosome/assembly annotation)

    Inputs:
        - assembly_id:      Sample identifier
        - assembly:         Unmasked genome assembly FASTA
                            (NOT soft-masked — centromeric satellites are tandem
                            repeats that RepeatMasker would mask, and centroAnno
                            needs to detect them via k-mer frequency analysis)

    Outputs:
        - centroanno_dir:   Full centroAnno output directory
        - monomer_csv:      Concatenated monomer decomposition CSV (all scaffolds)
        - hor_csv:          Concatenated HOR decomposition CSV (all scaffolds)
        - monomer_fa:       Concatenated monomer template FASTA (all scaffolds)
        - hor_fa:           Concatenated HOR template FASTA (all scaffolds)

    Tools: centroAnno (conda: qjh96::centroanno)

    Reference:
        Qi et al. (2025) De novo annotation of centromere with centroAnno.
        bioRxiv 2025.02.19.639205
========================================================================================
*/

process CENTROMERE_ANNOTATION {
    tag "${assembly_id}"
    label 'centromere'
    publishDir "${params.outdir}/centromere", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('_centroanno') ? null : fn }

    input:
    tuple val(assembly_id), path(assembly)

    output:
    tuple val(assembly_id), path("${assembly_id}_centroanno"),
        emit: centroanno_dir
    tuple val(assembly_id), path("${assembly_id}.centromere_monomers.csv"),
        emit: monomer_csv
    tuple val(assembly_id), path("${assembly_id}.centromere_hors.csv"),
        emit: hor_csv
    tuple val(assembly_id), path("${assembly_id}.monomerTemplates.fa"),
        emit: monomer_fa
    tuple val(assembly_id), path("${assembly_id}.HORs.fa"),
        emit: hor_fa

    script:
    """
    set -euo pipefail

    # ── Run centroAnno in assembly annotation mode ────────────────────────────
    centroAnno \\
        ${assembly} \\
        -o ${assembly_id}_centroanno \\
        -x anno-asm \\
        -t ${task.cpus}

    # ── Concatenate per-scaffold outputs into single files ────────────────────
    #   centroAnno creates per-sequence subdirectories under the output dir.
    #   Each subdir contains *decomposedResult.csv, *horDecomposedResult.csv,
    #   *monomerTemplates.fa, and *HORs.fa.

    # Monomer decomposition: seqname, template, start, end, identity, size
    > ${assembly_id}.centromere_monomers.csv
    find ${assembly_id}_centroanno -name '*decomposedResult.csv' \
        ! -name '*horDecomposedResult.csv' \
        -type f | sort | while read -r f; do
        cat "\$f" >> ${assembly_id}.centromere_monomers.csv
    done

    # HOR decomposition: seqname, template, start, end, identity, size, score
    > ${assembly_id}.centromere_hors.csv
    find ${assembly_id}_centroanno -name '*horDecomposedResult.csv' \
        -type f | sort | while read -r f; do
        cat "\$f" >> ${assembly_id}.centromere_hors.csv
    done

    # Monomer consensus templates
    > ${assembly_id}.monomerTemplates.fa
    find ${assembly_id}_centroanno -name '*monomerTemplates.fa' \
        -type f | sort | while read -r f; do
        cat "\$f" >> ${assembly_id}.monomerTemplates.fa
    done

    # HOR consensus sequences
    > ${assembly_id}.HORs.fa
    find ${assembly_id}_centroanno -name '*HORs.fa' \
        -type f | sort | while read -r f; do
        cat "\$f" >> ${assembly_id}.HORs.fa
    done
    """

    stub:
    """
    mkdir -p ${assembly_id}_centroanno
    touch ${assembly_id}.centromere_monomers.csv
    touch ${assembly_id}.centromere_hors.csv
    touch ${assembly_id}.monomerTemplates.fa
    touch ${assembly_id}.HORs.fa
    """
}

/*
========================================================================================
    CENTROANNO_TO_GFF3 MODULE
========================================================================================
    Convert centroAnno monomer and HOR decomposition CSVs to GFF3 format.

    centroAnno decomposedResult.csv columns (comma-separated):
        0: seqname      — scaffold/chromosome name
        1: template     — monomer template ID
        2: start        — 0-based start position
        3: end          — 0-based end position (exclusive)
        4: identity     — RepBlock-MonoTemp identity (0.0–1.0)
        5: size         — repeat unit size (bp)

    centroAnno horDecomposedResult.csv columns (comma-separated):
        0: seqname      — scaffold/chromosome name
        1: template     — HOR template ID
        2: start        — 0-based start position
        3: end          — 0-based end position (exclusive)
        4: identity     — RepBlock-MonoTemp identity (0.0–1.0)
        5: size         — HOR size (bp)
        6: score        — quality indicator (low/high)

    GFF3 output uses:
        - "satellite_DNA" (SO:0000005) for monomer repeat units
        - "repeat_region" (SO:0000657) for HOR annotations
        Source: centroAnno

    Inputs:
        - assembly_id:  Sample identifier
        - monomer_csv:  Concatenated monomer decomposition CSV
        - hor_csv:      Concatenated HOR decomposition CSV

    Outputs:
        - centromere_gff3:  ${assembly_id}.centromere.gff3

    Tools: Python 3 (stdlib only)
========================================================================================
*/

process CENTROMERE_TO_GFF3 {
    tag "${assembly_id}"
    label 'mergeFiles'
    publishDir "${params.outdir}/centromere", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(monomer_csv), path(hor_csv)

    output:
    tuple val(assembly_id), path("${assembly_id}.centromere.gff3"),
        emit: centromere_gff3

    script:
    """
    #!/usr/bin/env python3
    \"\"\"
    Convert centroAnno decomposition CSVs to GFF3 format.
    \"\"\"
    import csv
    import sys
    from collections import defaultdict

    gff_lines = []
    mono_count = 0
    hor_count = 0

    # ── Parse monomer decomposition ──────────────────────────────────────────
    #   Columns: seqname, template, start, end, identity, size
    #   Coordinates are 0-based; GFF3 requires 1-based inclusive.
    with open("${monomer_csv}") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split(",")
            if len(cols) < 6:
                continue

            seqname  = cols[0].strip()
            template = cols[1].strip()
            start_0  = int(cols[2].strip())
            end_0    = int(cols[3].strip())
            identity = cols[4].strip()
            size     = cols[5].strip()

            # Convert to 1-based inclusive coordinates
            start_1 = start_0 + 1
            end_1   = end_0  # end is exclusive in 0-based = inclusive in 1-based

            # Validate
            if end_1 < start_1:
                continue

            mono_count += 1

            # Score: use identity (0-1) scaled to 0-1000
            try:
                score = "{:.1f}".format(float(identity) * 1000)
            except ValueError:
                score = "."

            attrs = (
                f"ID=centromere_monomer_{mono_count};"
                f"Name={template};"
                f"repeat_unit_size={size};"
                f"identity={identity}"
            )

            gff_lines.append((
                seqname, "centroAnno", "satellite_DNA",
                start_1, end_1, score, ".", ".", attrs
            ))

    # ── Parse HOR decomposition ──────────────────────────────────────────────
    #   Columns: seqname, template, start, end, identity, size, score_class
    with open("${hor_csv}") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split(",")
            if len(cols) < 7:
                continue

            seqname     = cols[0].strip()
            template    = cols[1].strip()
            start_0     = int(cols[2].strip())
            end_0       = int(cols[3].strip())
            identity    = cols[4].strip()
            size        = cols[5].strip()
            score_class = cols[6].strip()

            start_1 = start_0 + 1
            end_1   = end_0

            if end_1 < start_1:
                continue

            hor_count += 1

            try:
                score = "{:.1f}".format(float(identity) * 1000)
            except ValueError:
                score = "."

            attrs = (
                f"ID=centromere_hor_{hor_count};"
                f"Name={template};"
                f"hor_unit_size={size};"
                f"identity={identity};"
                f"score_class={score_class}"
            )

            gff_lines.append((
                seqname, "centroAnno", "repeat_region",
                start_1, end_1, score, ".", ".", attrs
            ))

    # ── Sort by seqname, start, end ──────────────────────────────────────────
    gff_lines.sort(key=lambda x: (x[0], x[3], x[4]))

    # ── Write GFF3 ───────────────────────────────────────────────────────────
    with open("${assembly_id}.centromere.gff3", "w") as out:
        out.write("##gff-version 3\\n")
        for rec in gff_lines:
            out.write("\\t".join(str(x) for x in rec) + "\\n")

    print(
        f"centroAnno→GFF3: {mono_count} monomers, {hor_count} HORs "
        f"across {len(set(r[0] for r in gff_lines))} scaffolds",
        file=sys.stderr
    )
    """

    stub:
    """
    touch ${assembly_id}.centromere.gff3
    """
}