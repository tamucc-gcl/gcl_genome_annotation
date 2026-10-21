#!/usr/bin/env nextflow
// Development entry: completed selection, accepted before-QC and native post-filter QC.
nextflow.enable.dsl = 2
include { FINAL_STRUCTURAL_ANNOTATION } from '../../workflows/finalStructuralAnnotation.nf'

workflow {
    def root = file(params.structural_inputs_root ?: params.outdir).toString()
    def preflight = new groovy.json.JsonSlurperClassic().parse(new File("${root}/preflight/resolved_inputs.json"))
    if (preflight.status!='ready') error 'Final structural QC requires accepted input metadata'
    if (params.run_finalization) {
        def metadata = Channel.fromList(preflight.assemblies).map { row ->
            def aid = row.assembly_id
            if (!aid) error 'Missing structural QC assembly ID'
            tuple(aid,[run_busco:row.settings.run_busco,run_omark:row.settings.run_omark,
                busco_lineage:row.settings.busco_lineage,taxid:row.settings.taxid],
                file("${root}/repeats/${aid}.masked.fasta",checkIfExists:true),
                file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json",checkIfExists:true))
        }
        def selected = Channel.fromList(preflight.assemblies).map { row ->
            tuple(row.assembly_id,file("${root}/structural_annotation/selection/${row.assembly_id}",checkIfExists:true))
        }
        def assessments = Channel.fromList(preflight.assemblies).map { row ->
            def folder = "${root}/structural_annotation/assessment/${row.assembly_id}"
            tuple(row.assembly_id,file("${folder}/transcript_assessment.tsv",checkIfExists:true),file("${folder}/assessment_status.json",checkIfExists:true))
        }
        def beforeQc = params.run_baseline_qc ? Channel.fromList(preflight.assemblies).flatMap { row ->
            def aid = row.assembly_id
            def folder = "${root}/qc/final/${aid}"
            def files = ["${folder}/selection/selection.json","${folder}/agat/${aid}_annotation_stats"]
            if (row.settings.run_busco) files += "${folder}/busco/${aid}_busco.txt"
            if (row.settings.run_omark) files += "${folder}/omark/${aid}.summary.txt"
            files.collect { tuple(aid,file(it,checkIfExists:true)) }
        } : Channel.empty()
        FINAL_STRUCTURAL_ANNOTATION(metadata,selected,assessments,beforeQc)
    }
}
