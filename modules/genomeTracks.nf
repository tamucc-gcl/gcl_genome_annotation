process VALIDATE_GENOME_TRACK {
    tag "${aid}/${track}"
    label 'genome_track_report'
    publishDir "${params.outdir}/genome_tracks/${aid}/raw", mode: params.publish_dir_mode, overwrite: true
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

process HIGH_CONFIDENCE_TRNA {
    tag "${aid}"
    label 'trna_confidence'
    publishDir "${params.outdir}/genome_tracks/${aid}/raw", mode: params.publish_dir_mode, overwrite: true,
        saveAs: { name -> name=='confidence' ? 'trna_confidence' : name }
    input:
    tuple val(aid), path(table), path(structure)
    output:
    tuple val(aid), path('confidence'), emit: classifications
    script:
    """
    set -euo pipefail
    EukHighConfidenceFilter --result '${table}' --ss '${structure}' --output confidence --prefix '${aid}' \
      --cmscore1 ${params.trna_filter_cmscore} --ssscore1 ${params.trna_filter_ssscore} \
      --isoscore1 ${params.trna_filter_isoscore} --isoscore2 ${params.trna_filter_tertiary_start} \
      --isomaxscore2 ${params.trna_filter_tertiary_max}
    """
}

process FILTER_GENOME_TRACK {
    tag "${aid}/${track}"
    label 'genome_track_report'
    publishDir "${params.outdir}/genome_tracks/${aid}", mode: params.publish_dir_mode, overwrite: true
    input:
    tuple val(aid), val(track), path(source, stageAs:'incoming/raw.gff3'), path(raw_status, stageAs:'incoming/native.status.json'), path(identity), val(policy), path(aux, stageAs:'aux/*')
    path helpers
    output:
    tuple val(aid), val(track), path("${aid}.${track}.gff3"), path("${track}.status.json"), path("${track}_filter"), emit: filtered
    script:
    def encoded=groovy.json.JsonOutput.toJson(policy).getBytes('UTF-8').encodeBase64().toString()
    // A one-file path input is a TaskPath (itself iterable by path component).
    def auxiliaryFiles=aux instanceof Collection ? aux : [aux]
    def extra=auxiliaryFiles.collect { "'${it}'" }.join(' ')
    """
    python filter_genome_tracks.py --aid '${aid}' --track '${track}' --source '${source}' --identity '${identity}' \
      --native-status '${raw_status}' --policy '${encoded}' --aux ${extra} --outdir .
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
    def statusFiles=statuses instanceof Collection ? statuses : [statuses]
    def arguments = statusFiles.collect { "'${it}'" }.join(' ')
    """
    python report_genome_tracks.py --aid '${aid}' --settings '${encoded}' --identity '${identity}' --statuses ${arguments} --outdir .
    """
}
