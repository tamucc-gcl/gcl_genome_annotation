// Parsers are retained in the helper; accepted scientific metrics are joined unchanged.
process MERGE_FUNCTIONAL_ANNOTATIONS {
    tag "${assembly_id}"
    label 'functional_python'
    publishDir "${params.outdir}/functional_annotation/${assembly_id}", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), val(settings), path(assessment), path(eligibility),
          path(interproscan_tsv), path(eggnog_tsv), path(diamond_tsv), path(kofamscan_tsv), path(deepsig_gff3), path(omark_ump)
    path helpers
    output:
    tuple val(assembly_id), path("${assembly_id}.functional_annotations.tsv"), emit: merged_tsv
    tuple val(assembly_id), path('functional_merge_status.json'), emit: status
    script:
    def omarkArg = omark_ump ? "--omark '${omark_ump}'" : ''
    """
    python merge_functional_annotations.py --assembly-id '${assembly_id}' --assessment '${assessment}' \
      --eligibility '${eligibility}' --settings '${settings}' --interproscan-tsv '${interproscan_tsv}' \
      --eggnog-tsv '${eggnog_tsv}' --diamond-tsv '${diamond_tsv}' --kofamscan-tsv '${kofamscan_tsv}' \
      --deepsig-gff3 '${deepsig_gff3}' ${omarkArg} --outdir .
    """
}
