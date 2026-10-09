process AUDIT_REPEAT_LIBRARY {
    tag "${assembly_id}"
    publishDir "${params.outdir}/reference_repeats/assemblies", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), path(library)
    path helper
    output:
    tuple val(assembly_id), path("${assembly_id}.library.json"), emit: status
    script:
    """
    python ${helper} library --source '${library}' --assembly ${assembly_id} --output ${assembly_id}.library.json
    """
}

process VALIDATE_REFERENCE_REPEATS {
    tag "${assembly_id}"
    publishDir "${params.outdir}/reference_repeats/assemblies", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), path(reference_json), path(masked), path(annotations), path(library)
    val settings
    path helper
    output:
    tuple val(assembly_id), path("${assembly_id}.repeats.json"), emit: status
    script:
    """
    python ${helper} validate --reference-json '${reference_json}' --masked '${masked}' \
        --annotations '${annotations}' --library '${library}' --assembly ${assembly_id} --output ${assembly_id}.repeats.json
    python - <<'PY'
    import base64, json
    from pathlib import Path
    path = Path('${assembly_id}.repeats.json')
    data = json.loads(path.read_text())
    data['settings'] = json.loads(base64.b64decode('${settings}'))
    path.write_text(json.dumps(data, indent=2) + '\\n')
    PY
    """
}

process COMPARE_REPEAT_MASKS {
    tag "${assembly_id}"
    publishDir "${params.outdir}/reference_repeats/benchmarks", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), path(masked), path(annotations), path(unsplit_masked), path(unsplit_annotations)
    path helper
    output:
    path "${assembly_id}.benchmark.json"
    script:
    """
    python ${helper} compare --assembly ${assembly_id} --masked '${masked}' --annotations '${annotations}' \
        --files '${unsplit_masked}' '${unsplit_annotations}' --output ${assembly_id}.benchmark.json
    """
}

process REPORT_REFERENCE_REPEATS {
    publishDir "${params.outdir}/reference_repeats", mode: 'copy', overwrite: true
    input:
    path preflight
    path references, stageAs: 'references/reference??.json'
    path libraries, stageAs: 'libraries/library??.json'
    path repeats, stageAs: 'repeats/repeats??.json'
    path benchmarks, stageAs: 'benchmarks/benchmark??.json'
    path helper
    output:
    path 'reference_repeat_status.json', emit: status
    path 'README.md', emit: report
    script:
    """
    python ${helper} --preflight '${preflight}/resolved_inputs.json' --references references --libraries libraries --repeats repeats --benchmarks benchmarks --outdir .
    """
}

process ASSERT_REFERENCE_REPEATS {
    input:
    path status
    path benchmarks, stageAs: 'benchmarks/*'
    output:
    path 'reference_repeats_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    if json.loads(Path('${status}').read_text())['status'] != 'ready':
        raise SystemExit('Reference/repeat preparation incomplete; review reference_repeats/README.md')
    for path in Path('benchmarks').glob('*.json'):
        if json.loads(path.read_text())['status'] != 'ready':
            raise SystemExit('Whole-contig masking differs from the unsplit benchmark: ' + str(path))
    Path('reference_repeats_passed.txt').write_text('Reference identity and repeat outputs validated.\\n')
    PY
    """
}
