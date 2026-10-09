#!/usr/bin/env nextflow
// Foundation entry only. Does not launch annotation tools.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def request = [
        launch_dir: launchDir.toString(),
        sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(),
        taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }
    ]
    def encoded = groovy.json.JsonOutput.toJson(request).getBytes('UTF-8').encodeBase64().toString()
    RESOLVE_ANNOTATION_INPUTS(encoded,
        file("${projectDir}/py_scripts/resolve_annotation_inputs.py", checkIfExists: true),
        file("${projectDir}/assets/taxonomy_rules.json", checkIfExists: true))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
}
