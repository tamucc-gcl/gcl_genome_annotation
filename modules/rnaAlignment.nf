process PLAN_RNA_ALIGNMENT {
    cache false
    input:
    path preflight
    path passed
    path prepared
    val locations
    path helper
    output:
    path 'alignment_plan.json', emit: plan
    script:
    """
    python ${helper} plan --preflight '${preflight}/resolved_inputs.json' --prepared '${prepared}' --request '${locations}' --output alignment_plan.json
    """
}

process INDEX_RNA_REFERENCE {
    tag "${reference_id}"
    label 'rna_index'
    publishDir "${params.outdir}/rna_alignment/indices", mode: params.publish_dir_mode
    input:
    tuple val(reference_id), path(reference), path(identity)
    path helpers
    output:
    tuple val(reference_id), path("${reference_id}"), emit: index
    script:
    """
    set -euo pipefail
    mkdir ${reference_id}
    python rna_alignment.py reference --reference '${reference}' --identity '${identity}' --output '${reference_id}/reference.json'
    hisat2-build -p ${task.cpus} '${reference}' '${reference_id}/genome' > '${reference_id}/build.log' 2>&1
    hisat2-build --version > '${reference_id}/hisat2_version.txt'
    """
}

process ALIGN_RNA_READSET {
    tag "${aid}/${rid}"
    label 'rna_align'
    publishDir "${params.outdir}/rna_alignment/readsets/${aid}", mode: params.publish_dir_mode
    input:
    tuple val(aid), val(rid), val(spec), path(index), path(reads, stageAs: 'reads/mate??/*')
    path helpers
    output:
    tuple val(aid), val(rid), path("${rid}"), emit: alignment
    script:
    def options = new groovy.json.JsonSlurper().parseText(new String(spec.decodeBase64(), 'UTF-8'))
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def sampleArg = quote("SM:${options.sample_id ?: options.library_id}")
    def paths = reads instanceof Collection ? reads : [reads]
    def readArgs = options.layout == 'paired' ? "-1 ${quote(paths[0])} -2 ${quote(paths[1])}" : "-U ${quote(paths[0])}"
    def readCheck = groovy.json.JsonOutput.toJson([spec: options, paths: paths.collect { it.toString() }]).getBytes('UTF-8').encodeBase64().toString()
    def alignThreads = Math.min(params.rna_alignment_threads as int, Math.max(1, (task.cpus as int) - 2))
    def spliceReuseArg = params.rna_no_temp_splicesite ? '--no-temp-splicesite' : ''
    """
    set -euo pipefail
    python rna_alignment.py reads --request '${readCheck}'
    mkdir ${rid}
    # Strand-neutral mapping: splice-derived XS tags remain independent of a supplied strand setting.
    hisat2 -x '${index}/genome' ${readArgs} --dta ${spliceReuseArg} --seed ${params.rna_alignment_seed} \\
        --min-intronlen ${params.rna_min_intron} --max-intronlen ${params.rna_max_intron} \\
        --rg-id '${rid}' --rg ${sampleArg} --rg 'LB:${options.library_id}' \\
        -p ${alignThreads} 2> >(tee '${rid}/hisat2.log' >&2) | \\
        samtools sort -@ 1 -m 2G -o '${rid}/aligned.bam' -
    samtools index -c -@ ${alignThreads} '${rid}/aligned.bam'
    samtools quickcheck -v '${rid}/aligned.bam'
    samtools flagstat -@ ${alignThreads} '${rid}/aligned.bam' > '${rid}/flagstat.txt'
    samtools stats -@ ${alignThreads} '${rid}/aligned.bam' > '${rid}/samtools_stats.txt'
    hisat2 --version > '${rid}/hisat2_version.txt'
    samtools --version > '${rid}/samtools_version.txt'
    """
}

process QC_RNA_ALIGNMENT {
    tag "${aid}/${rid}"
    label 'rna_alignment_qc'
    publishDir "${params.outdir}/rna_alignment/qc/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(rid), val(spec), path(alignment), path(index)
    path helper
    output:
    tuple val(aid), val(rid), path("${rid}"), emit: qc
    script:
    """
    python ${helper} stats --request '${spec}' --bam '${alignment}/aligned.bam' --identity '${index}/reference.json' --outdir '${rid}'
    """
}

process RESOLVE_RNA_STRAND {
    tag "${aid}/${rid}"
    label 'rna_strand'
    publishDir "${params.outdir}/rna_alignment/strand/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(rid), val(policy), path(qc)
    path helper
    output:
    tuple val(aid), val(rid), path("${rid}.strand.json"), emit: strand
    script:
    """
    python ${helper} strand --request '${policy}' --identity '${qc}/alignment.json' --output '${rid}.strand.json'
    """
}

process REPORT_RNA_ALIGNMENT {
    label 'rna_strand'
    publishDir "${params.outdir}/rna_alignment", mode: 'copy', overwrite: true
    input:
    path plan
    path statistics, stageAs: 'qc/qc??.json'
    path strands, stageAs: 'strands/strand??.json'
    val bams
    path helper
    output:
    path 'rna_alignment_status.json', emit: status
    path 'README.md', emit: report
    script:
    """
    python ${helper} report --plan '${plan}' --qc qc --strands strands --request '${bams}' --outdir .
    """
}

process ASSERT_RNA_ALIGNMENT {
    label 'rna_strand'
    input:
    path status
    output:
    path 'rna_alignment_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    if json.loads(Path('${status}').read_text())['status'] != 'ready':
        raise SystemExit('RNA alignment incomplete; review rna_alignment/README.md')
    Path('rna_alignment_passed.txt').write_text('RNA alignment and readset BAMs validated.\\n')
    PY
    """
}
