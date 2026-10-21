#!/usr/bin/env nextflow
// Development entry: QC of accepted predictions, without repeating earlier stages.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from '../../modules/annotationInputs.nf'
include { ANNOTATION_BASELINE } from '../../workflows/annotationBaseline.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def root = file(params.baseline_inputs_root ?: params.outdir).toString()
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${params.pipeline_root}/py_scripts/resolve_annotation_inputs.py"),
        file("${params.pipeline_root}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    def accepted = RESOLVE_ANNOTATION_INPUTS.out.bundle.combine(ASSERT_ANNOTATION_INPUTS.out)
    def models = accepted.flatMap { folder, passed ->
        def preflight = new groovy.json.JsonSlurper().parse(new File("${folder}/resolved_inputs.json"))
        def report = new groovy.json.JsonSlurper().parse(new File("${root}/structural_annotation/prediction_status.json"))
        if (report.status != 'ready') error 'Baseline QC requires completed BRAKER predictions, not software checks'
        preflight.assemblies.collect { row ->
            def aid = row.assembly_id
            def prediction = report.assemblies.find { it.assembly_id == aid }
            if (!prediction || prediction.status != 'ready' || prediction.mode != row.prediction_mode)
                error "Missing/mismatched accepted prediction for ${aid}"
            def settings = [run_busco: row.settings.run_busco, run_omark: row.settings.run_omark,
                busco_lineage: row.settings.busco_lineage, taxid: row.settings.taxid]
            tuple(aid, settings, file("${root}/repeats/${aid}.masked.fasta", checkIfExists: true),
                file("${root}/structural_annotation/braker/${aid}/${aid}.gff3", checkIfExists: true),
                file("${root}/structural_annotation/braker/${aid}/${aid}_proteins.fa", checkIfExists: true),
                file("${root}/structural_annotation/braker_status/${aid}.prediction.json", checkIfExists: true))
        }
    }
    ANNOTATION_BASELINE(models, ['annotation_baseline.py', 'prepare_protein_evidence.py', 'reference_repeats.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
    })
}
