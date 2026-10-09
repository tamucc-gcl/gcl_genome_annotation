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
    // Keep real reads as path inputs. Local validation and passthrough stages emit
    // metadata only; Nextflow stages the original reads directly for each consumer.
    def rawReads = raw.ready.join(readsets).map { rid, bundle, row ->
        def metadata = new groovy.json.JsonSlurper().parse(new File("${bundle}/rna.json"))
        def paths = row.raw_request.run_accession ? metadata.files.collect { file("${bundle}/${it}") } :
            [row.fastq_1, row.fastq_2].findAll { it }.collect { file(it) }
        tuple(rid, file("${bundle}/rna.json"), paths)
    }
    def qcInput = { stage, rid, manifest, reads ->
        def suffixes = reads.collect { it.name.endsWith('.gz') ? '.fastq.gz' : '.fastq' }
        tuple(rid, encode([readset_id: rid, stage: stage]), suffixes + ['.fastq'], reads[0], reads.size() == 2 ? reads[1] : [])
    }
    QC_RNA_RAW(rawReads.map { rid, manifest, reads -> qcInput('raw', rid, manifest, reads) }, helper)
    def trimming = rawReads.join(readsets).map { rid, manifest, reads, row ->
        tuple(rid, encode([readset_id: rid] + row.trim_request), manifest, reads)
    }
    TRIM_RNA_READSET(trimming, helper)
    def clean = splitStatus(TRIM_RNA_READSET.out.bundle)
    def resolveReads = { rid, bundle, priorManifest, priorReads ->
        def metadata = new groovy.json.JsonSlurper().parse(new File("${bundle}/rna.json"))
        tuple(rid, file("${bundle}/rna.json"), metadata.passthrough ? priorReads :
            (metadata.status == 'ready' ? metadata.files.collect { file("${bundle}/${it}") } : []))
    }
    def cleanReads = clean.ready.join(rawReads).map(resolveReads)
    QC_RNA_PROCESSED(cleanReads.map { rid, manifest, reads -> qcInput('processed', rid, manifest, reads) }, helper)
    def selecting = cleanReads.join(readsets).map { rid, manifest, reads, row ->
        tuple(rid, encode([readset_id: rid] + row.subset_request), manifest, reads)
    }
    SUBSET_RNA_READSET(selecting, helper)
    def selectedReads = SUBSET_RNA_READSET.out.bundle.join(cleanReads).map(resolveReads)
    def failedReads = raw.failed.mix(clean.failed).map { rid, bundle -> tuple(rid, file("${bundle}/rna.json"), []) }
    def finalStates = selectedReads.mix(failedReads).toList()
    def finalResults = finalStates.map { rows -> rows.collect { it[1] } }
    def readPaths = finalStates.map { rows -> encode(rows.collectEntries { rid, manifest, reads ->
        [rid, reads.collect { it.toString() }]
    }) }
    // Stage individual ZIPs so MultiQC does not need to traverse linked directories.
    def qcResults = QC_RNA_RAW.out.qc.mix(QC_RNA_PROCESSED.out.qc).flatMap { rid, bundle ->
        def reports = file("${bundle}/*_fastqc.zip", checkIfExists: true)
        reports instanceof Collection ? reports : [reports]
    }.toList()
    def fastpReports = clean.ready.flatMap { rid, bundle ->
        def report = file("${bundle}/${rid}_fastp.json")
        report.exists() ? [report] : []
    }.toList()
    REPORT_RNA_PREPARATION(PLAN_RNA_PREPARATION.out.bundle, finalResults, readPaths, qcResults, fastpReports,
        file("${projectDir}/py_scripts/report_rna_preparation.py"))
    ASSERT_RNA_PREPARATION(REPORT_RNA_PREPARATION.out.status)
}
