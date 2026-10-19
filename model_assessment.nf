#!/usr/bin/env nextflow
// Assessment-only entry: accepted model/evidence outputs; no prediction or filtering.
nextflow.enable.dsl = 2
include { MODEL_ASSESSMENT } from './workflows/modelAssessment.nf'

workflow {
    def root = file(params.assessment_inputs_root ?: params.outdir).toString()
    def predictions = new groovy.json.JsonSlurperClassic().parse(new File("${root}/structural_annotation/prediction_status.json"))
    if (predictions.status != 'ready') error 'Model assessment requires accepted prediction metadata'
    if (params.run_model_assessment) {
        def rows = predictions.assemblies
        def annotations = Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            if (!aid || row.status != 'ready') error "Invalid assessment assembly: ${aid}"
            tuple(aid,file("${root}/repeats/${aid}.masked.fasta",checkIfExists:true),
                file("${root}/structural_annotation/final/${aid}/${aid}.gff3",checkIfExists:true),
                file("${root}/structural_annotation/final/${aid}/${aid}_proteins.fa",checkIfExists:true),
                file("${root}/structural_annotation/final/${aid}/${aid}_cds.fa",checkIfExists:true),
                file("${root}/structural_annotation/final/${aid}/annotation_status.json",checkIfExists:true),
                file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json",checkIfExists:true))
        }
        def repeats = params.run_te_flagging ? Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            def library = params.assessment_repeat_library ?: params.repeat_library
            if (!library) {
                def candidates = ["${root}/repeats/denovo/${aid}_reclassified_repeats.fa", "${root}/repeats/denovo/${aid}_curatedRepeats.fasta"]
                library = candidates.find { file(it).exists() }
            }
            if (!library) error "Missing pre-filter TE library for ${aid}; provide --assessment_repeat_library"
            tuple(aid,file("${root}/repeats/${aid}.repeat.gff3",checkIfExists:true),file(library,checkIfExists:true))
        } : Channel.empty()
        def aeds = params.run_aed ? Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            tuple(aid,file("${root}/structural_annotation/aed/${aid}/aed_scores.tsv",checkIfExists:true),
                file("${root}/structural_annotation/aed/${aid}/aed_status.json",checkIfExists:true))
        } : Channel.empty()
        def junctions = params.run_junction_qc ? Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            tuple(aid,file("${root}/qc/junctions/${aid}/transcript_junction_support.tsv",checkIfExists:true),
                file("${root}/qc/junctions/${aid}/junction_status.json",checkIfExists:true))
        } : Channel.empty()
        MODEL_ASSESSMENT(annotations,repeats,aeds,junctions)
    }
}
