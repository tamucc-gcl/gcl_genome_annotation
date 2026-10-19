nextflow.enable.dsl = 2
include { ANNOTATION_BASELINE as FINAL_ANNOTATION_QC } from './workflows/annotationBaseline.nf'
include { FINAL_STRUCTURAL_ANNOTATION } from './workflows/finalStructuralAnnotation.nf'

workflow {
    def root = file(params.structural_inputs_root).toString()
    def preflight = new groovy.json.JsonSlurperClassic().parse(new File("${root}/preflight/resolved_inputs.json"))
    if (params.run_finalization) {
        def rows = preflight.assemblies
        def metadata = Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            tuple(aid,row.settings,file("${root}/repeats/${aid}.masked.fasta"),file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json"))
        }
        FINAL_ANNOTATION_QC(metadata.map { aid,settings,genome,identity ->
            def folder = "${root}/structural_annotation/final/${aid}"
            tuple(aid,settings,genome,file("${folder}/${aid}.gff3"),file("${folder}/${aid}_proteins.fa"),file("${folder}/annotation_status.json"))
        },['annotation_baseline.py','prepare_protein_evidence.py','reference_repeats.py'].collect { file("${projectDir}/py_scripts/${it}") })
        def beforeQc = FINAL_ANNOTATION_QC.out.agat.mix(FINAL_ANNOTATION_QC.out.busco).mix(FINAL_ANNOTATION_QC.out.omark)
            .mix(FINAL_ANNOTATION_QC.out.selections.map { aid,folder -> tuple(aid,file("${folder}/selection.json")) })
        def selected = Channel.fromList(rows).map { row -> tuple(row.assembly_id,file("${root}/structural_annotation/selection/${row.assembly_id}")) }
        def assessments = Channel.fromList(rows).map { row ->
            def folder = "${root}/structural_annotation/assessment/${row.assembly_id}"
            tuple(row.assembly_id,file("${folder}/transcript_assessment.tsv"),file("${folder}/assessment_status.json"))
        }
        FINAL_STRUCTURAL_ANNOTATION(metadata,selected,assessments,beforeQc)
    }
}
