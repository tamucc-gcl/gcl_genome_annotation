#!/usr/bin/env nextflow
// Selection-only development entry; accepted assessment and models, no scientific searches.
nextflow.enable.dsl = 2
include { MODEL_SELECTION } from './workflows/modelSelection.nf'

workflow {
    def root = file(params.selection_inputs_root ?: params.outdir).toString()
    def predictions = new groovy.json.JsonSlurperClassic().parse(new File("${root}/structural_annotation/prediction_status.json"))
    if (predictions.status!='ready') error 'Model selection requires accepted prediction metadata'
    if (params.run_model_selection) {
        def rows = predictions.assemblies
        def annotations = Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            if (!aid || row.status!='ready') error "Invalid selection assembly: ${aid}"
            def folder = "${root}/structural_annotation/final/${aid}"
            tuple(aid,file("${folder}/${aid}.gff3",checkIfExists:true),file("${folder}/${aid}_proteins.fa",checkIfExists:true),
                file("${folder}/${aid}_cds.fa",checkIfExists:true),file("${folder}/annotation_status.json",checkIfExists:true),
                file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json",checkIfExists:true))
        }
        def assessments = Channel.fromList(rows).map { row ->
            def folder = "${root}/structural_annotation/assessment/${row.assembly_id}"
            tuple(row.assembly_id,file("${folder}/transcript_assessment.tsv",checkIfExists:true),
                file("${folder}/assessment_status.json",checkIfExists:true))
        }
        MODEL_SELECTION(annotations,assessments)
    }
}
