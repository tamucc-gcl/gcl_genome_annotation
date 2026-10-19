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
    output:
    path 'README.md',emit: report
    script:
    def expected = groovy.json.JsonOutput.toJson(ids).bytes.encodeBase64().toString()
    """
    python - <<'PY'
    import base64, json
    from pathlib import Path
    rows=[json.loads(p.read_text()) for p in Path('reports').glob('*/structural_qc_status.json')]
    expected=json.loads(base64.b64decode('${expected}'))
    if sorted(r['assembly_id'] for r in rows)!=sorted(expected): raise ValueError('Missing/duplicate assembly structural reports')
    lines=['# Structural annotation QC', '', '| Assembly | Before genes | After genes | Before transcripts | After transcripts |', '|---|---:|---:|---:|---:|']
    for row in sorted(rows,key=lambda r:r['assembly_id']):
        aid=row['assembly_id'];before,after=row['model_sets']
        lines.append(f"| [{aid}]({aid}/report.html) | {before['genes']} | {after['genes']} | {before['transcripts']} | {after['transcripts']} |")
    Path('README.md').write_text('\n'.join(lines)+'\n')
    PY
    """
}
