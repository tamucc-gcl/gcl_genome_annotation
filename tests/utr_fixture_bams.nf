nextflow.enable.dsl = 2

process CREATE_UTR_FIXTURE_BAMS {
    tag "${aid}"
    label 'rna_align'
    publishDir "${params.outdir}/rna_alignment", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(sam)
    output:
    path "readsets/${aid}/R1"
    script:
    """
    set -euo pipefail
    mkdir -p readsets/${aid}/R1
    samtools sort -o readsets/${aid}/R1/aligned.bam '${sam}'
    samtools index -c readsets/${aid}/R1/aligned.bam
    """
}

workflow {
    CREATE_UTR_FIXTURE_BAMS(Channel.fromPath("${params.fixture_sams}/*.sam").map { path -> tuple(path.baseName,path) })
}
