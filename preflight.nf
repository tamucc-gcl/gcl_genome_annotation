#!/usr/bin/env nextflow
// Foundation entry only. Does not include main.nf or launch annotation tools.
nextflow.enable.dsl = 2

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

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def request = [
        launch_dir: launchDir.toString(),
        sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(),
        taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }
    ]
    def encoded = groovy.json.JsonOutput.toJson(request).getBytes('UTF-8').encodeBase64().toString()
    RESOLVE_ANNOTATION_INPUTS(encoded,
        file("${projectDir}/py_scripts/resolve_annotation_inputs.py", checkIfExists: true),
        file("${projectDir}/assets/taxonomy_rules.json", checkIfExists: true))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
}
