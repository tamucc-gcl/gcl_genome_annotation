process VALIDATE_GENOME_TRACK {
    tag "${aid}/${track}"
    label 'genome_track_report'
    publishDir "${params.outdir}/genome_tracks/${aid}", mode: params.publish_dir_mode
    input:
    tuple val(aid), val(track), path(source), path(identity)
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
    publishDir "${params.outdir}/genome_tracks/${aid}", mode: params.publish_dir_mode
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
