// Original BRAKER3 training/prediction extended with matched, separate RNA BAMs.
process PLAN_BRAKER_INPUTS {
    cache false
    input:
    path preflight
    path repeats
    path proteins
    path rna
    path helpers
    output:
    path 'braker_plan.json', emit: plan
    script:
    """
    python braker_prediction.py plan --preflight '${preflight}/resolved_inputs.json' --repeats '${repeats}' \\
        --protein-status '${proteins}' --rna-status '${rna}' --output braker_plan.json
    """
}

process VALIDATE_BRAKER_INPUTS {
    tag "${aid}"
    label 'braker_inputs'
    publishDir "${params.outdir}/structural_annotation/braker_inputs", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(spec), path(genome), path(proteins), path(bams, stageAs: 'rna/readset??.bam'), path(indices, stageAs: 'rna/readset??.bam.csi'), path(identity)
    path helpers
    output:
    tuple val(aid), path("${aid}.inputs.json"), emit: status
    script:
    def bamList = bams instanceof Collection ? bams : [bams]
    def request = groovy.json.JsonOutput.toJson([genome: genome.toString(), proteins: proteins.toString(),
        bams: bamList.collect { it.toString() }, indices: (indices instanceof Collection ? indices : [indices]).collect { it.toString() }]).getBytes('UTF-8').encodeBase64().toString()
    """
    python braker_prediction.py validate-inputs --spec '${spec}' --request '${request}' --identity '${identity}' --output '${aid}.inputs.json'
    """
}

process CHECK_BRAKER_SOFTWARE {
    tag "${aid}"
    label 'braker_runtime'
    publishDir "${params.outdir}/structural_annotation/braker_runtime", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(spec), path(genome), path(proteins), path(bams, stageAs: 'rna/readset??.bam'), path(indices, stageAs: 'rna/readset??.bam.csi'), path(identity), path(validated)
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: runtime
    script:
    def bamList = bams instanceof Collection ? bams : [bams]
    def request = groovy.json.JsonOutput.toJson([genome: genome.toString(), proteins: proteins.toString(),
        bams: bamList.collect { it.toString() }, threads: 1, rounds: params.braker_rounds as int,
        workdir: 'software_check', augustus: 'augustus_config', ab_initio: false,
        check: true, container: params.braker_container]).getBytes('UTF-8').encodeBase64().toString()
    """
    set -euo pipefail
    mkdir augustus_config
    cp -a '${params.braker_augustus_config}/.' augustus_config/
    chmod -R u+rwX augustus_config
    python3 braker_prediction.py run --spec '${spec}' --request '${request}' --output '${aid}'
    """
}

process BRAKER {
    tag "${aid}"
    label 'braker'
    publishDir "${params.outdir}/structural_annotation/braker", mode: params.publish_dir_mode
    input:
    tuple val(aid), val(spec), path(genome), path(proteins), path(bams, stageAs: 'rna/readset??.bam'), path(indices, stageAs: 'rna/readset??.bam.csi'), path(identity), path(runtime)
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: prediction
    tuple val(aid), path("${aid}/${aid}.gff3"), emit: braker_gff3
    tuple val(aid), path("${aid}/${aid}_proteins.fa"), emit: braker_proteins
    tuple val(aid), path("${aid}/${aid}_cds.fa"), emit: braker_cds
    script:
    def threads = Math.min(params.braker_threads as int, task.cpus as int)
    def bamList = bams instanceof Collection ? bams : [bams]
    def request = groovy.json.JsonOutput.toJson([genome: genome.toString(), proteins: proteins.toString(),
        bams: bamList.collect { it.toString() }, threads: threads, rounds: params.braker_rounds as int,
        workdir: 'braker', augustus: 'augustus_config', ab_initio: params.braker_ab_initio,
        check: false, container: params.braker_container]).getBytes('UTF-8').encodeBase64().toString()
    def species = new groovy.json.JsonSlurper().parseText(new String(spec.decodeBase64(), 'UTF-8')).species
    """
    set -euo pipefail
    mkdir augustus_config
    cp -a '${params.braker_augustus_config}/.' augustus_config/
    chmod -R u+rwX augustus_config
    python3 braker_prediction.py run --spec '${spec}' --request '${request}' --output '${aid}'
    mv braker/braker.gff3 '${aid}/${aid}.gff3'
    mv braker/braker.aa '${aid}/${aid}_proteins.fa'
    mv braker/braker.codingseq '${aid}/${aid}_cds.fa'
    cp -a 'augustus_config/species/${species}' '${aid}/augustus_training'
    cp braker/braker.log '${aid}/braker.log'
    mv braker '${aid}/braker_run'
    """
}

process VALIDATE_BRAKER_OUTPUTS {
    tag "${aid}"
    label 'braker_report'
    publishDir "${params.outdir}/structural_annotation/braker_status", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(spec), path(prediction), path(identity)
    path helpers
    output:
    tuple val(aid), path("${aid}.prediction.json"), emit: status
    script:
    """
    python braker_prediction.py validate-outputs --spec '${spec}' --identity '${identity}' \\
        --gff '${prediction}/${aid}.gff3' --proteins '${prediction}/${aid}_proteins.fa' \\
        --cds '${prediction}/${aid}_cds.fa' --output '${aid}.prediction.json'
    """
}

process REPORT_BRAKER {
    label 'braker_report'
    publishDir "${params.outdir}/structural_annotation", mode: 'copy', overwrite: true
    input:
    path plan
    path predictions, stageAs: 'predictions/status??.json'
    path runtimes, stageAs: 'runtimes/runtime??/inputs.json'
    val check_only
    val rna_qc_link
    path helper
    output:
    path 'prediction_status.json', emit: status
    path 'README.md', emit: report
    script:
    """
    python ${helper} --plan '${plan}' --predictions predictions --runtimes runtimes --check-only '${check_only}' --rna-qc-link '${rna_qc_link}' --outdir .
    """
}

process ASSERT_BRAKER {
    label 'braker_report'
    input:
    path status
    output:
    path 'prediction_passed.txt'
    script:
    """
    python - <<'PY'
    import json
    from pathlib import Path
    state = json.loads(Path('${status}').read_text())
    if state['status'] not in ('ready', 'software_checked'):
        raise SystemExit('Prediction handoff incomplete; review structural_annotation/README.md')
    Path('prediction_passed.txt').write_text(state['status'] + '\\n')
    PY
    """
}
