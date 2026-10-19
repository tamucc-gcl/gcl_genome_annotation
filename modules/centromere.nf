/*
========================================================================================
    centroAnno MODULE
========================================================================================
    De novo centromere and tandem repeat annotation using centroAnno

    centroAnno identifies centromeric satellite arrays and their higher-order
    repeat (HOR) structure directly from genome assemblies without requiring
    prior knowledge of monomer templates.

    Mode: genome assembly (-G true), centroAnno 1.0.2 (whole-chromosome/assembly annotation)

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
    // Allocation can reserve the node without multiplying decomposition memory by 64 threads.
    def threads = Math.min(task.cpus as int, params.centromere_threads as int)
    if (threads < 1) error 'centromere_cpus and centromere_threads must be positive'
    """
    set -euo pipefail

    # ── Run centroAnno in assembly annotation mode ────────────────────────────
    centroAnno \\
        ${assembly} \\
        -o ${assembly_id}_centroanno \\
        -G true \\
        -t ${threads}

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

// Native CSV conversion retains coordinates and rejects unexpected schemas.
process CENTROMERE_TO_GFF3 {
    tag "${assembly_id}"
    label 'genome_track_report'
    publishDir "${params.outdir}/centromere", mode: params.publish_dir_mode
    input:
    tuple val(assembly_id), path(monomer_csv), path(hor_csv)
    path helpers
    output:
    tuple val(assembly_id), path("${assembly_id}.centromere.gff3"), emit: centromere_gff3
    script:
    """
    python genome_tracks.py centromeres --aid '${assembly_id}' --monomers '${monomer_csv}' --hors '${hor_csv}' --output '${assembly_id}.centromere.gff3'
    """
    stub:
    """
    touch ${assembly_id}.centromere.gff3
    """
}
