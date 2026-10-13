process PREPARE_ANNOTATION_BASELINE {
    tag "${aid}"
    label 'baseline_python'
    publishDir "${params.outdir}/qc/baseline/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(settings), path(genome), path(gff), path(proteins), path(validation)
    path helpers
    output:
    tuple val(aid), path('selection'), emit: selection
    script:
    def request = groovy.json.JsonOutput.toJson([assembly_id: aid, settings: settings,
        genome: genome.toString(), gff: gff.toString(), proteins: proteins.toString(), validation: validation.toString()])
        .getBytes('UTF-8').encodeBase64().toString()
    """
    python annotation_baseline.py prepare --request '${request}' --outdir selection
    """
}

process REPORT_ANNOTATION_BASELINE {
    label 'baseline_python'
    publishDir "${params.outdir}/qc/baseline", mode: 'copy', overwrite: true
    input:
    val ids
    path selections, stageAs: 'selection/*/selection.json'
    path agat, stageAs: 'agat/*/result.txt'
    val busco_ids
    path busco, stageAs: 'busco/*/result.txt'
    val omark_ids
    path omark, stageAs: 'omark/*/result.txt'
    path helpers
    output:
    path 'baseline_status.json', emit: status
    path 'README.md', emit: report
    script:
    def list = { value -> value instanceof Collection ? value : [value] }
    // stageAs's numbered wildcard keeps identical basenames distinct. Explicit
    // assembly IDs are carried separately and supplied as named input mappings.
    def request = groovy.json.JsonOutput.toJson([
        selections: list(selections).collect { it.toString() },
        agat: list(agat).withIndex().collect { p, i -> [aid: ids[i], path: p.toString()] },
        busco: list(busco).withIndex().collect { p, i -> [aid: busco_ids[i], path: p.toString()] },
        omark: list(omark).withIndex().collect { p, i -> [aid: omark_ids[i], path: p.toString()] }])
        .getBytes('UTF-8').encodeBase64().toString()
    """
    python annotation_baseline.py report --request '${request}' --outdir .
    """
}
