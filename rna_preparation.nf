#!/usr/bin/env nextflow
// Input preflight + assembly-independent RNA preparation; stops before alignment/repeats.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'
include { PLAN_RNA_PREPARATION; ACQUIRE_RNA_READSET; QC_RNA_RAW; TRIM_RNA_READSET; QC_RNA_PROCESSED; SUBSET_RNA_READSET; REPORT_RNA_PREPARATION; ASSERT_RNA_PREPARATION } from './modules/rnaPreparation.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    def helper = file("${projectDir}/py_scripts/prepare_rna.py", checkIfExists: true)
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${projectDir}/py_scripts/resolve_annotation_inputs.py"),
        file("${projectDir}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    PLAN_RNA_PREPARATION(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out,
        encode([quality: params.rna_trim_quality, window: params.rna_trim_window,
                min_length: params.rna_min_read_length, max_n: params.rna_max_n_bases]), helper)

    def readsets = PLAN_RNA_PREPARATION.out.bundle.flatMap { folder ->
        def plan = new groovy.json.JsonSlurper().parse(new File("${folder}/plan.json"))
        plan.readsets.collect { row -> tuple(row.readset_id, row) }
    }
    def acquisition = readsets.map { rid, row ->
        tuple(rid, encode(row.raw_request), [row.fastq_1, row.fastq_2].findAll { it }.collect { file(it, checkIfExists: true) })
    }
    def splitStatus = { channel -> channel.branch {
        ready: new groovy.json.JsonSlurper().parse(new File("${it[1]}/rna.json")).status == 'ready'
        failed: true
    } }
    ACQUIRE_RNA_READSET(acquisition, helper)
    def raw = splitStatus(ACQUIRE_RNA_READSET.out.bundle)
    QC_RNA_RAW(raw.ready.map { rid, bundle -> tuple(rid, encode([readset_id: rid, stage: 'raw']), bundle) }, helper)
    def trimming = raw.ready.join(readsets).map { rid, bundle, row ->
        tuple(rid, encode([readset_id: rid] + row.trim_request), bundle)
    }
    TRIM_RNA_READSET(trimming, helper)
    def clean = splitStatus(TRIM_RNA_READSET.out.bundle)
    QC_RNA_PROCESSED(clean.ready.map { rid, bundle -> tuple(rid, encode([readset_id: rid, stage: 'processed']), bundle) }, helper)
    def selecting = clean.ready.join(readsets).map { rid, bundle, row ->
        tuple(rid, encode([readset_id: rid] + row.subset_request), bundle)
    }
    SUBSET_RNA_READSET(selecting, helper)
    def finalResults = SUBSET_RNA_READSET.out.bundle.mix(raw.failed, clean.failed).map { rid, bundle -> bundle }.toList()
    // Stage individual ZIPs so MultiQC does not need to traverse linked directories.
    def qcResults = QC_RNA_RAW.out.qc.mix(QC_RNA_PROCESSED.out.qc).flatMap { rid, bundle ->
        def reports = file("${bundle}/*_fastqc.zip", checkIfExists: true)
        reports instanceof Collection ? reports : [reports]
    }.toList()
    def fastpReports = clean.ready.flatMap { rid, bundle ->
        def report = file("${bundle}/${rid}_fastp.json")
        report.exists() ? [report] : []
    }.toList()
    REPORT_RNA_PREPARATION(PLAN_RNA_PREPARATION.out.bundle, finalResults, qcResults, fastpReports,
        file("${projectDir}/py_scripts/report_rna_preparation.py"))
    ASSERT_RNA_PREPARATION(REPORT_RNA_PREPARATION.out.status)
}
