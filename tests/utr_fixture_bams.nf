nextflow.enable.dsl = 2

process CREATE_UTR_FIXTURE_BAMS {
    tag "${aid}"
    label 'rna_align'
    publishDir "${params.outdir}/rna_alignment", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(sam)
    output:
    path "readsets/${aid}"
    script:
    """
    set -euo pipefail
    for rid in R1 R2; do
      mkdir -p readsets/${aid}/\$rid
      samtools sort -o readsets/${aid}/\$rid/aligned.bam '${sam}'
      samtools index -c readsets/${aid}/\$rid/aligned.bam
    done
    """
}

workflow {
    CREATE_UTR_FIXTURE_BAMS(Channel.fromPath("${params.fixture_sams}/*.sam").map { path -> tuple(path.baseName,path) })
}
