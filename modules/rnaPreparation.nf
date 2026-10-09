process PLAN_RNA_PREPARATION {
    cache false
    input:
    path preflight
    path passed
    val options
    path helper
    output:
    path 'rna_plan', emit: bundle
    script:
    """
    python ${helper} plan --manifest '${preflight}/resolved_inputs.json' --request '${options}' --outdir rna_plan
    """
}

process ACQUIRE_RNA_READSET {
    tag "${rid}"
    publishDir "${params.outdir}/rna/raw", mode: 'symlink', overwrite: true
    input:
    tuple val(rid), val(request), path(reads, stageAs: 'incoming/*')
    path helper
    output:
    tuple val(rid), path("${rid}"), emit: bundle
    script:
    def paths = reads instanceof Collection ? reads : [reads]
    def inputs = groovy.json.JsonOutput.toJson(paths.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${helper} acquire --request '${request}' --inputs '${inputs}' --outdir ${rid} --threads ${task.cpus}
    """
}

process QC_RNA_RAW {
    tag "${rid}"
    publishDir "${params.outdir}/rna/qc/raw", mode: 'copy', overwrite: true
    input:
    tuple val(rid), val(request), path(bundle, stageAs: 'raw_bundle')
    path helper
    output:
    tuple val(rid), path("${rid}_raw"), emit: qc
    script:
    """
    python ${helper} qc --request '${request}' --bundle ${bundle} --outdir ${rid}_raw --threads ${task.cpus}
    """
}

process TRIM_RNA_READSET {
    tag "${rid}"
    publishDir "${params.outdir}/rna/processed", mode: 'symlink', overwrite: true
    input:
    tuple val(rid), val(request), path(bundle, stageAs: 'raw_bundle')
    path helper
    output:
    tuple val(rid), path("${rid}"), emit: bundle
    script:
    """
    python ${helper} trim --request '${request}' --bundle ${bundle} --outdir ${rid} --threads ${task.cpus}
    """
}

process QC_RNA_PROCESSED {
    tag "${rid}"
    publishDir "${params.outdir}/rna/qc/processed", mode: 'copy', overwrite: true
    input:
    tuple val(rid), val(request), path(bundle, stageAs: 'processed_bundle')
    path helper
    output:
    tuple val(rid), path("${rid}_processed"), emit: qc
    script:
    """
    python ${helper} qc --request '${request}' --bundle ${bundle} --outdir ${rid}_processed --threads ${task.cpus}
    """
}

process SUBSET_RNA_READSET {
    tag "${rid}"
    publishDir "${params.outdir}/rna/selected", mode: 'symlink', overwrite: true
    input:
    tuple val(rid), val(request), path(bundle, stageAs: 'processed_bundle')
    path helper
    output:
    tuple val(rid), path("${rid}"), emit: bundle
    script:
    """
    python ${helper} subset --request '${request}' --bundle ${bundle} --outdir ${rid}
    """
}

process REPORT_RNA_PREPARATION {
    publishDir "${params.outdir}/rna", mode: 'copy', overwrite: true
    input:
    path plan
    path results, stageAs: 'results/*'
    path qc, stageAs: 'qc/*'
    path fastp_reports, stageAs: 'fastp/*'
    path helper
    output:
    path 'rna_status.json', emit: status
    path 'README.md', emit: report
    path 'multiqc', emit: qc
    script:
    """
    python ${helper} --plan ${plan}/plan.json --results results --qc qc --fastp fastp --outdir .
    """
}

process ASSERT_RNA_PREPARATION {
    input:
    path status
    output:
    path 'rna_preparation_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    data = json.loads(Path('${status}').read_text())
    if data['status'] != 'ready':
        raise SystemExit('RNA preparation incomplete. Review rna/README.md and readset manifests.')
    Path('rna_preparation_passed.txt').write_text('RNA reads prepared; reference alignment/auto strandedness remain pending.\\n')
    PY
    """
}
