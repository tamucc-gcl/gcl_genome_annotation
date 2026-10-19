process PREPARE_MODEL_METRICS {
    tag "${aid}"
    label 'aed_python'
    input:
    tuple val(aid), path(genome), path(gff), path(proteins), path(cds), path(validation), path(identity)
    path helpers
    output:
    tuple val(aid), path('metrics/metrics.json'), path('metrics/regions.bed'), emit: metrics
    script:
    """
    python model_metrics.py --assembly-id '${aid}' --genome '${genome}' --gff '${gff}' --proteins '${proteins}' --cds '${cds}' --validation '${validation}' --identity '${identity}' --outdir metrics
    """
}

process REPORT_MODEL_ASSESSMENT {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/assessment/${aid}", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(metrics), path(aed_files,stageAs:'aed/*'), path(junction_files,stageAs:'junctions/*'), path(te_files,stageAs:'te/*')
    path helpers
    output:
    tuple val(aid), path('assessment_status.json'), emit: status
    tuple val(aid), path('transcript_assessment.tsv'), emit: transcripts
    tuple val(aid), path('gene_assessment.tsv'), emit: genes
    tuple val(aid), path('README.md'), emit: report
    script:
    def options = (aed_files ? '--aed-status aed/aed_status.json --aed-scores aed/aed_scores.tsv' : '') + ' ' +
        (junction_files ? '--junction-status junctions/junction_status.json --junction-scores junctions/transcript_junction_support.tsv' : '') + ' ' +
        (te_files ? '--te-status te/te_status.json --te-scores te/te_gene_report.tsv' : '')
    """
    python model_assessment.py --assembly-id '${aid}' --metrics '${metrics}' ${options} --outdir .
    """
}
