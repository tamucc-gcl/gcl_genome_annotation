// Unified Markdown results/methods report. No searches or model changes.
process ANNOTATION_REPORT {
    label 'annotation_report'
    publishDir "${params.outdir}", mode:'copy', overwrite:true, pattern:'annotation_report.md'
    publishDir "${params.outdir}/reports", mode:'copy', overwrite:true, pattern:'annotation_report_*'
    publishDir "${params.outdir}/reports", mode:'copy', overwrite:true, pattern:'annotation_metrics.tsv'
    input:
    path releases, stageAs:'release_inputs/*'
    val structural_ids
    path structural, stageAs:'structural_inputs/*'
    val stage_keys
    path stages, stageAs:'stage??/*'
    val configuration
    path helpers
    output:
    path 'annotation_report.md', emit: report
    path 'annotation_report_plots', emit: plots
    path 'annotation_report_data', emit: data
    path 'annotation_metrics.tsv', emit: metrics
    path 'annotation_report_status.json', emit: status
    script:
    def asList = { value -> value instanceof Collection ? value : [value] }
    def structuralFiles=asList(structural)
    def stageFiles=asList(stages)
    def request=[releases:asList(releases).collect { it.toString() },
        structural:structural_ids.withIndex().collect { aid,i -> [aid,structuralFiles[i].toString()] },
        stages:stage_keys.withIndex().collect { key,i -> [key[0],key[1],stageFiles[i].toString()] },configuration:configuration]
    def encoded=groovy.json.JsonOutput.toJson(request).bytes.encodeBase64().toString()
    """
    python generate_annotation_report.py --request '${encoded}' --outdir .
    """
}
