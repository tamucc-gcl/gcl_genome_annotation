// Shared preflight processes extracted unchanged from the cluster-validated entry.
process RESOLVE_ANNOTATION_INPUTS {
    cache false  // Also checks original shared-filesystem paths for presence on each invocation.
    // The declared output is already a directory named preflight.
    // Publish the bundle under outdir, yielding outdir/preflight exactly once.
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    input:
    val request_base64
    path resolver
    path rules
    output:
    path 'preflight', emit: bundle
    script:
    """
    python -c "import base64; open('request.json','wb').write(base64.b64decode('${request_base64}'))"
    python ${resolver} --request request.json --rules ${rules} --outdir preflight --report-only
    """
}

process ASSERT_ANNOTATION_INPUTS {
    input:
    path bundle
    output:
    path 'preflight_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    manifest = json.loads(Path('${bundle}/resolved_inputs.json').read_text())
    if manifest['status'] != 'ready':
        raise SystemExit('Annotation preflight incomplete. Review the published preflight README.md and resolved_inputs.json.')
    Path('preflight_passed.txt').write_text('Input configuration passed; biological/read-content checks remain pending.\\n')
    PY
    """
}

