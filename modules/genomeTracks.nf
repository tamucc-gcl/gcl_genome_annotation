process VALIDATE_GENOME_TRACK {
    tag "${aid}/${track}"
    label 'genome_track_report'
    publishDir "${params.outdir}/genome_tracks/${aid}", mode: params.publish_dir_mode, overwrite: true
    input:
    // Keep the staged native file distinct from the normalized published output.
    tuple val(aid), val(track), path(source, stageAs: 'incoming/native.gff3'), path(identity)
    path helpers
    output:
    tuple val(aid), val(track), path("${aid}.${track}.gff3"), path("${track}.status.json"), emit: validated
    script:
    """
    python genome_tracks.py normalize --aid '${aid}' --track '${track}' --source '${source}' --identity '${identity}' --output '${aid}.${track}.gff3'
    """
}

process REPORT_GENOME_TRACKS {
    tag "${aid}"
    label 'genome_track_report'
    // Resume must replace a report published by a different enabled/disabled policy.
    publishDir "${params.outdir}/genome_tracks/${aid}", mode: params.publish_dir_mode, overwrite: true
    input:
    tuple val(aid), val(settings), path(identity), path(statuses)
    path helpers
    output:
    tuple val(aid), path('genome_tracks_status.json'), emit: status
    path('README.md'), emit: readme
    script:
    def encoded = groovy.json.JsonOutput.toJson(settings).getBytes('UTF-8').encodeBase64().toString()
    def arguments = statuses.collect { "'${it}'" }.join(' ')
    """
    python report_genome_tracks.py --aid '${aid}' --settings '${encoded}' --identity '${identity}' --statuses ${arguments} --outdir .
    """
}
