#!/usr/bin/env nextflow
// Input preflight + normalized reference/repeats; stops before RNA alignment/prediction.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'
include { FILTER_ASSEMBLY } from './modules/filterAssembly.nf'
include { CURATE_DFAM } from './modules/curateDfam.nf'
include { REPEAT_MODEL } from './modules/repeatModel.nf'
include { CURATE_REPEATS } from './modules/curateRepeats.nf'
include { CLASSIFY_REPEATS } from './modules/classifyRepeats.nf'
include { RECLASSIFY_REPEATS } from './modules/reclassifyRepeats.nf'
include { EXTRACT_REPEAT_PEPS; INDEX_REPEAT_PEPS; FILTER_REPEAT_LIBRARY } from './modules/filterRepeatLibrary.nf'
include { COMBINE_REPEAT_LIBS } from './modules/combineRepeatLibs.nf'
include { REPEAT_MASK_PARALLEL; REPEAT_MASK_PARALLEL as UNSPLIT_MASK_BENCHMARK } from './workflows/repeatMaskParallel.nf'
include { AUDIT_REPEAT_LIBRARY; VALIDATE_REFERENCE_REPEATS; COMPARE_REPEAT_MASKS; REPORT_REFERENCE_REPEATS; ASSERT_REFERENCE_REPEATS } from './modules/referenceRepeats.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    if (params.n_rm_chunks < 1) error '--n_rm_chunks must be positive'
    if (params.repeat_seed < 0) error '--repeat_seed must be nonnegative'
    if (!(params.job_speed in [null, 's', 'q', 'qq'])) error '--job_speed must be s, q, qq or null'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def parse = { path -> new groovy.json.JsonSlurper().parse(new File(path.toString())) }
    def taxonKey = { taxa -> 'taxon_' + java.security.MessageDigest.getInstance('SHA-256')
        .digest(taxa.toString().getBytes('UTF-8')).encodeHex().toString().substring(0, 20) }
    def helper = file("${projectDir}/py_scripts/reference_repeats.py", checkIfExists: true)
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${projectDir}/py_scripts/resolve_annotation_inputs.py"),
        file("${projectDir}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    // Assertion dependency ensures no scientific work starts on an incomplete input contract.
    def assemblies = RESOLVE_ANNOTATION_INPUTS.out.bundle.combine(ASSERT_ANNOTATION_INPUTS.out).flatMap { folder, passed ->
        parse("${folder}/resolved_inputs.json").assemblies.collect { row ->
            if (params.repeat_mask_benchmark && row.assembly_id.endsWith('__benchmark'))
                error 'Assembly IDs ending __benchmark are reserved when the masking benchmark is enabled'
            tuple(row.assembly_id, row)
        }
    }
    FILTER_ASSEMBLY(assemblies.map { aid, row ->
        tuple(aid, file(row.assembly_fasta, checkIfExists: true), row.settings.min_contig_length)
    }, helper)
    def references = FILTER_ASSEMBLY.out.bundle.map { aid, folder ->
        tuple(aid, file("${folder}/${aid}.reference.json"))
    }
    def split = references.branch {
        ready: parse(it[1]).status == 'ready'
        failed: true
    }
    def prepared = split.ready.join(FILTER_ASSEMBLY.out.bundle).map { aid, state, folder ->
        tuple(aid, file("${folder}/${aid}.filtered.fasta"), state)
    }
    def librarySources
    if (params.repeat_library) {
        // An explicit complete library is a scientific input, also used by tiny masking fixtures.
        def supplied = file(params.repeat_library, checkIfExists: true)
        librarySources = prepared.map { aid, reference, state -> tuple(aid, supplied) }
    } else {
        if (!params.dfam_database) error 'Provide --dfam_database for the existing Dfam FamDB directory'
        if (!params.diamond_database && !params.repeat_filter_database) error 'Provide --diamond_database or --repeat_filter_database'
        def dfam = file("${params.dfam_database}/dfam_db", checkIfExists: true)
        def generalDb = file(params.repeat_filter_database ?: "${params.diamond_database}/swissprot.dmnd", checkIfExists: true)
        def requestedTaxa = prepared.join(assemblies).map { aid, reference, state, row ->
            tuple(taxonKey(row.settings.taxa_name), aid, row.settings.taxa_name)
        }
        CURATE_DFAM(requestedTaxa.map { key, aid, taxa -> tuple(key, taxa, dfam) }.unique { it[0] })
        def dfamByAssembly = requestedTaxa.combine(CURATE_DFAM.out.bundle, by: 0).map { key, aid, taxa, folder ->
            tuple(aid, file("${folder}/dfam_subset.lib"), file("${folder}/dfam_subset.source.json"))
        }
        REPEAT_MODEL(prepared.join(dfamByAssembly).map { aid, reference, state, library, source ->
            tuple(aid, reference, source)
        }, params.repeat_seed)
        CURATE_REPEATS(REPEAT_MODEL.out.families_fa, helper)
        CLASSIFY_REPEATS(CURATE_REPEATS.out.clustered_repeats)
        def classified
        if (params.repeat_reclassify) {
            RECLASSIFY_REPEATS(CLASSIFY_REPEATS.out.curated_repeats, file("${projectDir}/py_scripts/run_terrier.py"))
            classified = RECLASSIFY_REPEATS.out.reclassified_lib
        } else {
            classified = CLASSIFY_REPEATS.out.curated_repeats
        }
        EXTRACT_REPEAT_PEPS()
        INDEX_REPEAT_PEPS(EXTRACT_REPEAT_PEPS.out.repeat_peps.map { path ->
            tuple('repeat_peps_' + java.security.MessageDigest.getInstance('SHA-256').digest(path.bytes).encodeHex().toString(), path)
        })
        // No-input acquisition and its mapped index already emit reusable value channels.
        FILTER_REPEAT_LIBRARY(classified, generalDb, EXTRACT_REPEAT_PEPS.out.repeat_peps,
            INDEX_REPEAT_PEPS.out.index.map { key, db, source -> db })
        COMBINE_REPEAT_LIBS(FILTER_REPEAT_LIBRARY.out.filtered_lib.join(dfamByAssembly).map { aid, library, dfamLib, source ->
            tuple(aid, library, dfamLib)
        })
        librarySources = COMBINE_REPEAT_LIBS.out.combined_lib
    }
    AUDIT_REPEAT_LIBRARY(librarySources, helper)
    def acceptedLibraries = AUDIT_REPEAT_LIBRARY.out.status.filter { aid, status -> parse(status).status == 'ready' }
        .join(librarySources).map { aid, status, library -> tuple(aid, library) }
    def maskInputs = prepared.join(acceptedLibraries).map { aid, reference, state, library ->
        tuple(aid, reference, library)
    }
    REPEAT_MASK_PARALLEL(maskInputs, params.n_rm_chunks)
    def benchmarkChecks = Channel.value([])
    if (params.repeat_mask_benchmark) {
        UNSPLIT_MASK_BENCHMARK(maskInputs.map { aid, reference, library -> tuple(aid + '__benchmark', reference, library) }, 1)
        def benchmark = UNSPLIT_MASK_BENCHMARK.out.masked_fasta.join(UNSPLIT_MASK_BENCHMARK.out.repeat_out)
            .map { aid, masked, annotations -> tuple(aid.toString().replaceFirst(/__benchmark$/, ''), masked, annotations) }
        COMPARE_REPEAT_MASKS(REPEAT_MASK_PARALLEL.out.masked_fasta.join(REPEAT_MASK_PARALLEL.out.repeat_out)
            .join(benchmark), helper)
        benchmarkChecks = COMPARE_REPEAT_MASKS.out.toList()
    }
    def validations = REPEAT_MASK_PARALLEL.out.masked_fasta.join(REPEAT_MASK_PARALLEL.out.repeat_out)
        .join(split.ready).join(acceptedLibraries).map { aid, masked, annotations, state, library ->
            tuple(aid, state, masked, annotations, library)
        }
    VALIDATE_REFERENCE_REPEATS(validations, encode([
        tetools_version: params.tetools_version, dfam_version: params.dfam_version,
        repeat_seed: params.repeat_seed, additional_model_rounds: params.repModel_rounds,
        min_repeat_length: params.min_repeat_length, cluster_pct: params.cluster_pct,
        repeat_reclassify: params.repeat_reclassify, terrier_threshold: params.terrier_threshold,
        terrier_environment: 'bio-terrier=0.4.0; torch=2.6.0; lightning=2.6.1',
        filter_evalue: params.repeat_filter_evalue, filter_pident: params.repeat_filter_pident,
        library_mode: params.repeat_library ? 'supplied_complete_library' : 'denovo_plus_dfam',
        masking_policy: 'whole_contigs', requested_batches: params.n_rm_chunks, job_speed: params.job_speed
    ]), helper)
    REPORT_REFERENCE_REPEATS(RESOLVE_ANNOTATION_INPUTS.out.bundle, references.map { it[1] }.toList(),
        AUDIT_REPEAT_LIBRARY.out.status.map { it[1] }.toList(), VALIDATE_REFERENCE_REPEATS.out.status.map { it[1] }.toList(),
        benchmarkChecks,
        file("${projectDir}/py_scripts/report_reference_repeats.py"))
    ASSERT_REFERENCE_REPEATS(REPORT_REFERENCE_REPEATS.out.status, benchmarkChecks)
}
