// Structural QC adapted from the existing gene aggregation and report metrics.
process STRUCTURAL_QC {
    tag "${aid}"
    label 'structural_report'
    publishDir "${params.outdir}/qc/structural", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(identified,stageAs:'identified_source'), path(selection,stageAs:'selected_source'), path(before_table), path(before_status), path(before_qc,stageAs:'before_qc/*'), path(after_qc,stageAs:'after_qc/*')
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: report
    script:
    def beforeOption = before_qc ? "--before-qc ${before_qc}" : ''
    """
    python structural_report.py --assembly-id '${aid}' --identified '${identified}' --selection '${selection}' \
      --before-table '${before_table}' --before-status '${before_status}' ${beforeOption} --after-qc ${after_qc} --outdir '${aid}'
    """
}

process STRUCTURAL_QC_INDEX {
    label 'aed_python'
    publishDir "${params.outdir}/qc/structural", mode:'copy', overwrite:true
    input:
    val ids
    path reports,stageAs:'reports/*'
    path index_helper
    output:
    path 'README.md',emit: report
    script:
    def expected = groovy.json.JsonOutput.toJson(ids).bytes.encodeBase64().toString()
    """
    python structural_report_index.py --expected '${expected}' --reports ${reports} --out README.md
    """
}
