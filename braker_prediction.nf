#!/usr/bin/env nextflow
// Development entry: accepted real outputs -> BRAKER3, without repeats/remapping.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'
include { PLAN_BRAKER_INPUTS; VALIDATE_BRAKER_INPUTS; CHECK_BRAKER_SOFTWARE; BRAKER; VALIDATE_BRAKER_OUTPUTS; REPORT_BRAKER; ASSERT_BRAKER } from './modules/braker.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    if ((params.braker_threads as int) < 1 || (params.braker_threads as int) > 48 || (params.braker_cpus as int) < 1 || (params.braker_rounds as int) < 1)
        error 'BRAKER needs positive CPUs/rounds and 1..48 tool threads'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    def root = file(params.braker_inputs_root ?: params.outdir).toString()
    def helpers = ['braker_prediction.py', 'reference_repeats.py', 'prepare_protein_evidence.py'].collect {
        file("${projectDir}/py_scripts/${it}", checkIfExists: true)
    }
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${projectDir}/py_scripts/resolve_annotation_inputs.py"),
        file("${projectDir}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    def accepted = RESOLVE_ANNOTATION_INPUTS.out.bundle.combine(ASSERT_ANNOTATION_INPUTS.out).map { bundle, passed -> bundle }
    PLAN_BRAKER_INPUTS(accepted, file("${root}/reference_repeats/reference_repeat_status.json", checkIfExists: true),
        file("${root}/protein_evidence/protein_status.json", checkIfExists: true),
        file("${root}/rna_alignment/rna_alignment_status.json", checkIfExists: true), helpers)
    def jobs = PLAN_BRAKER_INPUTS.out.plan.flatMap { path -> new groovy.json.JsonSlurper().parse(path.toFile()).jobs }
    def inputs = jobs.map { row ->
        def aid = row.assembly_id
        tuple(aid, encode(row), file("${root}/repeats/${aid}.masked.fasta", checkIfExists: true),
            file("${root}/protein_evidence/assemblies/${aid}/proteins.fa", checkIfExists: true),
            row.readsets.collect { file("${root}/rna_alignment/${it.bam}", checkIfExists: true) },
            row.readsets.collect { file("${root}/rna_alignment/${it.bam_index}", checkIfExists: true) },
            file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json", checkIfExists: true))
    }
    VALIDATE_BRAKER_INPUTS(inputs, helpers)
    CHECK_BRAKER_SOFTWARE(inputs.join(VALIDATE_BRAKER_INPUTS.out.status), helpers)
    def results = Channel.value([])
    if (!params.braker_check_only) {
        BRAKER(inputs.join(CHECK_BRAKER_SOFTWARE.out.runtime), helpers)
        VALIDATE_BRAKER_OUTPUTS(inputs.join(BRAKER.out.prediction).map { aid, spec, genome, proteins, bams, indices, identity, prediction ->
            tuple(aid, spec, prediction, identity)
        }, helpers)
        results = VALIDATE_BRAKER_OUTPUTS.out.status.map { aid, path -> path }.toList()
    }
    REPORT_BRAKER(PLAN_BRAKER_INPUTS.out.plan, results,
        CHECK_BRAKER_SOFTWARE.out.runtime.map { aid, folder -> file("${folder}/inputs.json") }.toList(),
        params.braker_check_only,
        file("${params.outdir}/structural_annotation").toAbsolutePath().relativize(file("${root}/rna_alignment/multiqc/multiqc_report.html").toAbsolutePath()).toString(),
        file("${projectDir}/py_scripts/report_braker.py"))
    ASSERT_BRAKER(REPORT_BRAKER.out.status)
}
