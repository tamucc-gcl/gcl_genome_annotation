process PLAN_PROTEIN_EVIDENCE {
    cache false // Recheck source presence/identity and previously failed acquisitions.
    input:
    path preflight_bundle
    path preflight_passed
    val options_base64
    path helper
    path catalogue
    output:
    path 'plan', emit: bundle
    script:
    """
    python ${helper} plan --manifest '${preflight_bundle}/resolved_inputs.json' \
        --catalogue ${catalogue} --options-base64 '${options_base64}' --outdir plan
    """
}

process ACQUIRE_PROTEIN_SOURCE {
    tag "${source_id}"
    label 'protein_source'
    input:
    tuple val(source_id), val(spec_base64), path(local_inputs, stageAs: 'input/*')
    path helper
    output:
    tuple val(source_id), path("${source_id}"), emit: source
    script:
    def inputs = local_inputs instanceof Collection ? local_inputs : [local_inputs]
    def staged = groovy.json.JsonOutput.toJson(inputs.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python -c "import base64; open('source_request.json','wb').write(base64.b64decode('${spec_base64}'))"
    python ${helper} acquire --spec source_request.json --staged-base64 '${staged}' --outdir ${source_id}
    """
}

process PUBLISH_FAILED_PROTEIN_SOURCE {
    tag "${source_id}"
    publishDir "${params.outdir}/protein_sources", mode: 'symlink', overwrite: true
    input:
    tuple val(source_id), path(bundle, stageAs: 'incoming/*')
    output:
    tuple val(source_id), path("${source_id}"), emit: source
    script:
    """
    python - <<'PY'
    from pathlib import Path
    Path('${source_id}').symlink_to(Path('${bundle}').resolve(), target_is_directory=True)
    PY
    """
}

// Only successful acquisitions reach this process; failed bundles are never stored.
process STORE_PROTEIN_SOURCE {
    tag "${source_id}"
    storeDir "${store_directory}"
    publishDir "${params.outdir}/protein_sources", mode: 'symlink', overwrite: true
    input:
    tuple val(source_id), val(store_directory), path(bundle, stageAs: 'incoming/*')
    output:
    tuple val(source_id), path("${source_id}"), emit: source
    script:
    """
    python - <<'PY'
    import json, shutil
    from pathlib import Path
    source = Path('${bundle}')
    assert json.loads((source / 'source.json').read_text())['status'] == 'ready'
    shutil.copytree(source, '${source_id}')
    Path('${source_id}/STORE_READY.json').write_text('{"status":"ready","normalization_version":1}\\n')
    PY
    """
}

process COMBINE_PROTEIN_EVIDENCE {
    tag "${assembly_id}"
    label 'protein_combine'
    publishDir "${params.outdir}/protein_evidence/assemblies", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), val(assembly_base64), path(source_bundles, stageAs: 'sources/*')
    path helper
    output:
    tuple val(assembly_id), path("${assembly_id}"), emit: evidence
    script:
    """
    python -c "import base64; open('assembly_request.json','wb').write(base64.b64decode('${assembly_base64}'))"
    python ${helper} combine --assembly assembly_request.json --sources-directory sources --outdir ${assembly_id}
    """
}

process REPORT_PROTEIN_EVIDENCE {
    publishDir "${params.outdir}/protein_evidence", mode: 'copy', overwrite: true
    input:
    path evidence_bundles, stageAs: 'assemblies/*'
    path report_helper
    output:
    path 'protein_status.json', emit: status
    path 'README.md', emit: report
    script:
    def evidence = evidence_bundles instanceof Collection ? evidence_bundles : [evidence_bundles]
    def bundles = groovy.json.JsonOutput.toJson(evidence.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${report_helper} --bundles-base64 '${bundles}' \
        --status protein_status.json --report README.md
    """
}

process ASSERT_PROTEIN_EVIDENCE {
    input:
    path status
    output:
    path 'protein_preparation_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    if json.loads(Path('${status}').read_text())['status'] != 'ready':
        raise SystemExit('Protein preparation incomplete. Review protein_evidence/README.md and source error logs.')
    Path('protein_preparation_passed.txt').write_text('Requested protein sources validated and assembled.\\n')
    PY
    """
}
