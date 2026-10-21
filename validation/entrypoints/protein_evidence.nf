#!/usr/bin/env nextflow
// Input preflight + prediction protein preparation only. No RNA, repeats or BRAKER.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from '../../modules/annotationInputs.nf'
include { PLAN_PROTEIN_EVIDENCE; ACQUIRE_PROTEIN_SOURCE; PUBLISH_FAILED_PROTEIN_SOURCE; STORE_PROTEIN_SOURCE; COMBINE_PROTEIN_EVIDENCE; REPORT_PROTEIN_EVIDENCE; ASSERT_PROTEIN_EVIDENCE } from '../../modules/proteinEvidence.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def helper = file("${params.pipeline_root}/py_scripts/prepare_protein_evidence.py", checkIfExists: true)
    RESOLVE_ANNOTATION_INPUTS(encode(request),
        file("${params.pipeline_root}/py_scripts/resolve_annotation_inputs.py", checkIfExists: true),
        file("${params.pipeline_root}/assets/taxonomy_rules.json", checkIfExists: true))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    PLAN_PROTEIN_EVIDENCE(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out,
        encode([launch_dir: launchDir.toString(), orthodb_database: params.orthodb_database?.toString(),
                protein_source_store: params.protein_source_store?.toString(),
                published_sources: "${params.outdir}/protein_sources"]), helper,
        file("${params.pipeline_root}/assets/orthodb_sources.json", checkIfExists: true))

    def plan_sources = PLAN_PROTEIN_EVIDENCE.out.bundle.flatMap { directory ->
        def data = new groovy.json.JsonSlurper().parse(new File("${directory}/plan.json"))
        data.sources.collect { source ->
            // Create the shared store parent on the driver before acquisition.
            // A parent first created by a scratch worker can be temporarily
            // invisible to the driver's immediate storeDir output check.
            def root = source.kind == 'orthodb' ?
                "${params.orthodb_database}/validated/${source.release}/normalization-v1" :
                "${params.protein_source_store}/normalization-v1"
            java.nio.file.Files.createDirectories(file(root))
            tuple(source.source_id, encode(source),
                  source.local_path ? [file(source.local_path, checkIfExists: true)] : [])
        }
    }
    def plan_assemblies = PLAN_PROTEIN_EVIDENCE.out.bundle.flatMap { directory ->
        def data = new groovy.json.JsonSlurper().parse(new File("${directory}/plan.json"))
        data.assemblies.collect { assembly ->
            tuple(assembly.assembly_id, encode(assembly))
        }
    }
    def routed = plan_sources.branch {
        stored: new groovy.json.JsonSlurper().parseText(new String(it[1].decodeBase64(), 'UTF-8')).stored_path
        acquire: true
    }
    ACQUIRE_PROTEIN_SOURCE(routed.acquire, helper)
    def existing = routed.stored.map { sid, encoded, inputs ->
        def spec = new groovy.json.JsonSlurper().parseText(new String(encoded.decodeBase64(), 'UTF-8'))
        tuple(sid, file(spec.stored_path, checkIfExists: true))
    }
    def resolved = ACQUIRE_PROTEIN_SOURCE.out.source.mix(existing).branch {
        ready: new groovy.json.JsonSlurper().parse(new File("${it[1]}/source.json")).status == 'ready'
        failed: true
    }
    def stores = resolved.ready.map { sid, bundle ->
        // Store placement is transport metadata, independent of source identity.
        def plans = new groovy.json.JsonSlurper().parse(new File("${bundle}/request.json"))
        def root = plans.kind == 'orthodb' ?
            "${params.orthodb_database}/validated/${plans.release}/normalization-v1" :
            "${params.protein_source_store}/normalization-v1"
        tuple(sid, file(root).toString(), bundle)
    }
    STORE_PROTEIN_SOURCE(stores)
    PUBLISH_FAILED_PROTEIN_SOURCE(resolved.failed)
    def sources = STORE_PROTEIN_SOURCE.out.source.mix(PUBLISH_FAILED_PROTEIN_SOURCE.out.source).map { sid, bundle -> bundle }.toList()
    // Only each assembly's requested bundles enter its combination/cache identity.
    def combinations = plan_assemblies.combine(sources.map { bundles -> [bundles] }).map { aid, encoded, bundles ->
        def assembly = new groovy.json.JsonSlurper().parseText(new String(encoded.decodeBase64(), 'UTF-8'))
        tuple(aid, encoded, bundles.findAll { assembly.source_ids.contains(it.fileName.toString()) }.sort { it.fileName.toString() })
    }
    COMBINE_PROTEIN_EVIDENCE(combinations, helper)
    REPORT_PROTEIN_EVIDENCE(COMBINE_PROTEIN_EVIDENCE.out.evidence.map { aid, bundle -> bundle }.collect(),
        file("${params.pipeline_root}/py_scripts/report_protein_evidence.py", checkIfExists: true))
    ASSERT_PROTEIN_EVIDENCE(REPORT_PROTEIN_EVIDENCE.out.status)
}
