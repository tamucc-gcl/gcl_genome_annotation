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
    // Dynamic path names are resolved per task in Nextflow 23.10.1. A GString
    // inside a stageAs option map is evaluated while declaring the process.
    tuple val(rid), val(request), val(suffixes), path("${rid}_raw_R1${suffixes[0]}"), path("${rid}_raw_R2${suffixes.size() > 1 ? suffixes[1] : '.fastq'}")
    path helper
    output:
    tuple val(rid), path("${rid}_raw"), emit: qc
    script:
    def reads = (0..<suffixes.size()).collect { "${rid}_raw_R${it + 1}${suffixes[it]}" }
    def inputs = groovy.json.JsonOutput.toJson(reads.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${helper} qc --request '${request}' --inputs '${inputs}' --outdir ${rid}_raw --threads ${task.cpus}
    """
}

process TRIM_RNA_READSET {
    tag "${rid}"
    publishDir "${params.outdir}/rna/processed", mode: 'symlink', overwrite: true
    input:
    tuple val(rid), val(request), path(manifest, stageAs: 'raw.json'), path(reads, stageAs: 'reads/*')
    path helper
    output:
    tuple val(rid), path("${rid}"), emit: bundle
    script:
    def paths = reads instanceof Collection ? reads : [reads]
    def inputs = groovy.json.JsonOutput.toJson(paths.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${helper} trim --request '${request}' --bundle ${manifest} --inputs '${inputs}' --outdir ${rid} --threads ${task.cpus}
    """
}

process QC_RNA_PROCESSED {
    tag "${rid}"
    publishDir "${params.outdir}/rna/qc/processed", mode: 'copy', overwrite: true
    input:
    tuple val(rid), val(request), val(suffixes), path("${rid}_processed_R1${suffixes[0]}"), path("${rid}_processed_R2${suffixes.size() > 1 ? suffixes[1] : '.fastq'}")
    path helper
    output:
    tuple val(rid), path("${rid}_processed"), emit: qc
    script:
    def reads = (0..<suffixes.size()).collect { "${rid}_processed_R${it + 1}${suffixes[it]}" }
    def inputs = groovy.json.JsonOutput.toJson(reads.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${helper} qc --request '${request}' --inputs '${inputs}' --outdir ${rid}_processed --threads ${task.cpus}
    """
}

process SUBSET_RNA_READSET {
    tag "${rid}"
    publishDir "${params.outdir}/rna/selected", mode: 'symlink', overwrite: true
    input:
    tuple val(rid), val(request), path(manifest, stageAs: 'processed.json'), path(reads, stageAs: 'reads/*')
    path helper
    output:
    tuple val(rid), path("${rid}"), emit: bundle
    script:
    def paths = reads instanceof Collection ? reads : [reads]
    def inputs = groovy.json.JsonOutput.toJson(paths.collect { it.toString() }).getBytes('UTF-8').encodeBase64().toString()
    """
    python ${helper} subset --request '${request}' --bundle ${manifest} --inputs '${inputs}' --outdir ${rid}
    """
}

process REPORT_RNA_PREPARATION {
    publishDir "${params.outdir}/rna", mode: 'copy', overwrite: true
    input:
    path plan
    path results, stageAs: 'results/rna??.json'
    val read_paths
    path raw_qc, stageAs: 'raw_qc/*'
    path processed_qc, stageAs: 'processed_qc/*'
    path fastp_reports, stageAs: 'fastp/*'
    path helper
    output:
    path 'rna_status.json', emit: status
    path 'README.md', emit: report
    path 'multiqc', emit: qc
    script:
    """
    python ${helper} --plan ${plan}/plan.json --results results --read-paths '${read_paths}' --raw-qc raw_qc --processed-qc processed_qc --fastp fastp --outdir .
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
