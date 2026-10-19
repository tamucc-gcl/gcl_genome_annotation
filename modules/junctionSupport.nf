process RUN_PORTCULLIS {
    tag "${aid}/${rid}"
    label 'portcullis'
    publishDir "${params.outdir}/qc/junctions/${aid}/readsets", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(rid), val(request), val(strand), path(genome), path(identity), path(bam, stageAs: 'alignment.bam'), path(bai, stageAs: 'alignment.bam.csi')
    path helpers
    output:
    tuple val(aid), val(rid), path("${rid}"), emit: junctions
    script:
    """
    set -euo pipefail
    mkdir '${rid}'
    python portcullis_inputs.py check --genome '${genome}' --identity '${identity}' --bam '${bam}' --request '${request}' --out '${rid}/inputs.json'
    portcullis --version 2>&1 | tee '${rid}/version.txt'
    portcullis prep --use_csi -t ${task.cpus} -o prep '${genome}' '${bam}' 2>&1 | tee '${rid}/prep.log'
    portcullis junc --use_csi -t ${task.cpus} --strandedness '${strand}' --orientation UNKNOWN -o '${rid}/all' prep 2>&1 | tee '${rid}/analysis.log'
    portcullis filter -t ${task.cpus} --threshold ${params.portcullis_threshold} --save_bad --intron_gff -o '${rid}/filtered' prep '${rid}/all.junctions.tab' 2>&1 | tee '${rid}/filter.log'
    python portcullis_inputs.py finish --folder '${rid}' --threshold ${params.portcullis_threshold}
    """
}

process REPORT_JUNCTION_SUPPORT {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/qc/junctions/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(request), path(gff), path(validation), path(identity), path(bundles, stageAs: 'readsets/*')
    path helpers
    output:
    tuple val(aid), path('junction_status.json'), emit: status
    tuple val(aid), path('transcript_junction_support.tsv'), emit: transcripts
    tuple val(aid), path('junction_support.tsv'), emit: junctions
    tuple val(aid), path('README.md'), emit: report
    script:
    def paths = bundles instanceof Collection ? bundles : [bundles]
    def options = paths.collect { "--readset '${it}'" }.join(' ')
    """
    python junction_support.py --gff '${gff}' --validation '${validation}' --identity '${identity}' --request '${request}' ${options} --outdir .
    """
}
