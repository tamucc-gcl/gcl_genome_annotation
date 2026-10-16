process PREPARE_INGENANNOT_EVIDENCE {
    tag "${aid}"
    label 'ingenannot'
    input:
    tuple val(aid), path(gff), path(validation), path(current_status), path(rna, stageAs: 'rna/*'), path(protein, stageAs: 'protein/*')
    path helpers
    output:
    tuple val(aid), path('indexed'), emit: evidence
    script:
    def r = rna instanceof Collection ? rna : [rna]
    def p = protein instanceof Collection ? protein : [protein]
    def options = (r ? "--rna '${r[0]}'" : '') + ' ' + (p ? "--protein '${p[0]}'" : '')
    def rpath = r ? "'${r[0]}'" : '/dev/null'
    def ppath = p ? "'${p[0]}'" : '/dev/null'
    def prior = current_status instanceof Collection ? current_status : [current_status]
    def priorOption = prior ? "--current-status '${prior[0]}'" : ''
    """
    set -euo pipefail
    mkdir indexed
    python ingenannot_inputs.py --gff '${gff}' --validation '${validation}' \
      ${priorOption} ${options} --out indexed/inputs.json
    LC_ALL=C sort -k1,1 -k4,4n -k5,5n ${rpath} | bgzip -c > indexed/rna.gtf.gz
    tabix -p gff indexed/rna.gtf.gz
    LC_ALL=C sort -k1,1 -k4,4n -k5,5n ${ppath} | bgzip -c > indexed/protein.gff.gz
    tabix -p gff indexed/protein.gff.gz
    """
}

process RUN_INGENANNOT_AED {
    tag "${aid}"
    label 'ingenannot'
    publishDir "${params.outdir}/structural_annotation/${task.ext.aed_stage ?: 'ingenannot_benchmark'}/${aid}/native", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(gff), path(indexed)
    output:
    tuple val(aid), path('native.aed.gff3'), emit: annotated
    tuple val(aid), path('ingenannot_version.txt'), emit: version
    script:
    def cds = params.ingenannot_rna_cds_only ? '--aed_tr_cds_only' : ''
    """
    set -euo pipefail
    ingenannot -v 2 -p ${task.cpus} aed '${gff}' native.aed.gff3 braker3 \
      '${indexed}/rna.gtf.gz' '${indexed}/protein.gff.gz' \
      --evtr_source stringtie --evpr_source miniprot --evtrstranded --evprstranded \
      --penalty_overflow ${params.ingenannot_rna_penalty} ${cds}
    python -c "from importlib.metadata import version; print(version('ingenannot'))" > ingenannot_version.txt
    """
}

process REPORT_INGENANNOT_BENCHMARK {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/ingenannot_benchmark/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(gff), path(identity), path(indexed), path(current_scores), path(annotated)
    path helpers
    output:
    tuple val(aid), path('benchmark_scores.tsv'), emit: scores
    tuple val(aid), path('benchmark_status.json'), emit: status
    tuple val(aid), path('README.md'), emit: report
    script:
    """
    python ingenannot_benchmark.py report --original '${gff}' --annotated '${annotated}' \
      --identity '${identity}' --inputs '${indexed}/inputs.json' --current-scores '${current_scores}' \
      --cds-only '${params.ingenannot_rna_cds_only}' --penalty ${params.ingenannot_rna_penalty} --outdir .
    """
}
