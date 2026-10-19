process PREPARE_FUNCTIONAL_INPUTS {
    tag "${aid}"
    label 'functional_python'
    publishDir "${params.outdir}/functional_annotation/${aid}", mode: 'copy', overwrite: true, pattern: 'eligibility.tsv'
    input:
    tuple val(aid), path(proteins), path(gff), path(assessment), path(validation)
    path helpers
    output:
    tuple val(aid), path('eligible.fa'), path('eligibility.tsv'), path('inputs.json'), emit: prepared
    script:
    """
    python functional_inputs.py --assembly-id '${aid}' --proteins '${proteins}' --gff '${gff}' \
      --assessment '${assessment}' --validation '${validation}' --outdir .
    """
}

// A bound empty result for disabled tools or annotations with no eligible proteins.
// These are metadata/results placeholders, never links or copies of input data.
process EMPTY_FUNCTIONAL_RESULTS {
    tag "${aid}"
    label 'functional_python'
    input:
    val aid
    output:
    tuple val(aid), path('empty_results'), emit: results
    script:
    """
    mkdir empty_results
    touch empty_results/interproscan.tsv empty_results/eggnog.tsv empty_results/kofamscan.tsv empty_results/deepsig.gff3
    printf 'qseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\n' > empty_results/diamond.tsv
    """
}

process REPORT_FUNCTIONAL_ANNOTATION {
    tag "${aid}"
    label 'functional_python'
    publishDir "${params.outdir}/functional_annotation/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(original), path(annotated), path(table), path(status), path(eligibility)
    path helpers
    output:
    tuple val(aid), path('README.md'), emit: report
    tuple val(aid), path('functional_status.json'), emit: status
    script:
    """
    python functional_report.py --assembly-id '${aid}' --original '${original}' --annotated '${annotated}' \
      --table '${table}' --merge-status '${status}' --eligibility '${eligibility}' --outdir .
    """
}
