#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PREPARE_INGENANNOT_EVIDENCE; RUN_INGENANNOT_AED; REPORT_INGENANNOT_BENCHMARK } from './modules/ingenannotBenchmark.nf'

workflow {
    if (params.ingenannot_rna_penalty < 0 || params.ingenannot_rna_penalty > 1)
        error 'ingenannot_rna_penalty must be between 0 and 1'
    def root = file(params.aed_inputs_root ?: params.outdir).toString()
    def predictions = new groovy.json.JsonSlurper().parse(new File("${root}/structural_annotation/prediction_status.json"))
    if (predictions.status != 'ready') error 'Benchmark requires completed BRAKER predictions'
    def helpers = ['ingenannot_benchmark.py','evidence_aed.py','prepare_protein_evidence.py'].collect { file("${projectDir}/py_scripts/${it}") }
    def requests = Channel.fromList(predictions.assemblies).map { row ->
        def aid = row.assembly_id
        if (row.status != 'ready') error "Prediction not ready: ${aid}"
        def rna = row.mode in ['rna','rna_protein']
        def protein = row.mode in ['protein','rna_protein']
        tuple(aid, file("${root}/structural_annotation/braker/${aid}/${aid}.gff3", checkIfExists: true),
          file("${root}/structural_annotation/braker_status/${aid}.prediction.json", checkIfExists: true),
          file("${root}/structural_annotation/evidence_aed/${aid}/aed_status.json", checkIfExists: true),
          rna ? [file("${root}/structural_annotation/braker/${aid}/braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff", checkIfExists: true)] : [],
          protein ? [file("${root}/structural_annotation/evidence_aed/${aid}/protein_alignment/${aid}.miniprot.gff", checkIfExists: true)] : [])
    }
    PREPARE_INGENANNOT_EVIDENCE(requests, helpers)
    def ready = requests.map { aid,gff,validation,status,rna,protein -> tuple(aid,gff) }.join(PREPARE_INGENANNOT_EVIDENCE.out.evidence)
    RUN_INGENANNOT_AED(ready)
    REPORT_INGENANNOT_BENCHMARK(ready.join(RUN_INGENANNOT_AED.out.annotated).map { aid,gff,indexed,annotated ->
        tuple(aid,gff,file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json", checkIfExists: true),indexed,
          file("${root}/structural_annotation/evidence_aed/${aid}/aed_scores.tsv", checkIfExists: true),annotated)
    }, helpers)
}
