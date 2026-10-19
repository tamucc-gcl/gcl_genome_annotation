// Existing GFF decoration/product hierarchy, with optional other feature tracks.
process ANNOTATE_GFF3 {
    tag "${assembly_id}"
    label 'functional_python'
    publishDir "${params.outdir}/functional_annotation/${assembly_id}", mode: 'copy', overwrite: true
    input:
    tuple val(assembly_id), path(braker_gff3), path(trnascan_gff3), path(infernal_gff3),
          path(telomere_gff3), path(centromere_gff3), path(repeat_gff), path(functional_tsv)
    path helper
    output:
    tuple val(assembly_id), path("${assembly_id}.annotated.gff3"), emit: annotated_gff3
    tuple val(assembly_id), path("${assembly_id}.gene_product_summary.tsv"), emit: gene_summary
    script:
    def extraArgs = [[trnascan_gff3,'trnascan-gff3'],[infernal_gff3,'infernal-gff3'],
        [telomere_gff3,'telomere-gff3'],[centromere_gff3,'centromere-gff3'],[repeat_gff,'repeat-gff']]
        .findAll { it[0] }.collect { "--${it[1]} '${it[0]}'" }.join(' ')
    """
    python annotate_gff3.py --assembly-id '${assembly_id}' --braker-gff3 '${braker_gff3}' \
      --functional-tsv '${functional_tsv}' ${extraArgs}
    """
}
