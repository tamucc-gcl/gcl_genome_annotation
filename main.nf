#!/usr/bin/env nextflow
// Integrated rework: inputs -> proteins/RNA/reference/repeats -> alignment -> BRAKER3.
// Validated UTR refinement precedes final unfiltered QC and native AED assessment.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'
include { PLAN_PROTEIN_EVIDENCE; ACQUIRE_PROTEIN_SOURCE; PUBLISH_FAILED_PROTEIN_SOURCE; STORE_PROTEIN_SOURCE; COMBINE_PROTEIN_EVIDENCE; REPORT_PROTEIN_EVIDENCE; ASSERT_PROTEIN_EVIDENCE } from './modules/proteinEvidence.nf'
include { PLAN_RNA_PREPARATION; ACQUIRE_RNA_READSET; QC_RNA_RAW; TRIM_RNA_READSET; QC_RNA_PROCESSED; SUBSET_RNA_READSET; REPORT_RNA_PREPARATION; ASSERT_RNA_PREPARATION } from './modules/rnaPreparation.nf'
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
include { INDEX_RNA_REFERENCE; ALIGN_RNA_READSET; QC_RNA_ALIGNMENT; RESOLVE_RNA_STRAND; REPORT_RNA_ALIGNMENT; ASSERT_RNA_ALIGNMENT } from './modules/rnaAlignment.nf'
include { PLAN_BRAKER_INPUTS; VALIDATE_BRAKER_INPUTS; CHECK_BRAKER_SOFTWARE; BRAKER; VALIDATE_BRAKER_OUTPUTS; REPORT_BRAKER; ASSERT_BRAKER } from './modules/braker.nf'
include { PLAN_CONNECTED_RNA_ALIGNMENT } from './modules/connectedRnaAlignment.nf'
include { ANNOTATION_POSTPROCESSING } from './workflows/annotationPostprocessing.nf'
include { FUNCTIONAL_ANNOTATION } from './workflows/functionalAnnotation.nf'
include { GENOME_TRACKS } from './workflows/genomeTracks.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    if (params.run_functional_annotation && !params.run_finalization && !params.braker_check_only)
        error 'Functional annotation requires run_finalization=true (or disable run_functional_annotation)'
    if (params.n_rm_chunks < 1 || params.repeat_seed < 0) error 'Invalid repeat batches/seed'
    if (!(params.job_speed in [null, 's', 'q', 'qq'])) error '--job_speed must be s, q, qq or null'
    if ((params.rna_min_intron as int) < 1 || (params.rna_max_intron as int) < (params.rna_min_intron as int)) error 'Invalid RNA intron limits'
    if ((params.rna_strand_min_mapq as int) < 0 || (params.rna_strand_min_mapq as int) > 255) error 'Invalid RNA strand MAPQ'
    if ((params.rna_alignment_cpus as int) < 3 || (params.rna_alignment_threads as int) < 1) error 'Invalid RNA CPU/thread settings'
    if ((params.braker_threads as int) < 1 || (params.braker_threads as int) > 48 || (params.braker_cpus as int) < 1 || (params.braker_rounds as int) < 1)
        error 'BRAKER needs positive CPUs/rounds and 1..48 tool threads'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def parse = { path -> new groovy.json.JsonSlurper().parse(new File(path.toString())) }
    def taxonKey = { taxa -> 'taxon_' + java.security.MessageDigest.getInstance('SHA-256').digest(taxa.toString().getBytes('UTF-8')).encodeHex().toString().substring(0, 20) }
    def proteinHelper = file("${projectDir}/py_scripts/prepare_protein_evidence.py", checkIfExists: true)
    def rnaHelper = file("${projectDir}/py_scripts/prepare_rna.py", checkIfExists: true)
    def repeatHelper = file("${projectDir}/py_scripts/reference_repeats.py", checkIfExists: true)
    def alignmentHelper = file("${projectDir}/py_scripts/rna_alignment.py", checkIfExists: true)
    def helpers = [alignmentHelper, repeatHelper]
    def readHelpers = [alignmentHelper, rnaHelper]
    def brakerHelpers = ['braker_prediction.py', 'reference_repeats.py', 'prepare_protein_evidence.py'].collect {
        file("${projectDir}/py_scripts/${it}", checkIfExists: true)
    }
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${projectDir}/py_scripts/resolve_annotation_inputs.py"),
        file("${projectDir}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)

    // Protein evidence: shared database stores and project-specific spike-ins.
    PLAN_PROTEIN_EVIDENCE(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out,
        encode([launch_dir: launchDir.toString(), orthodb_database: params.orthodb_database?.toString(),
                protein_source_store: params.protein_source_store?.toString(),
                published_sources: "${params.outdir}/protein_sources"]), proteinHelper,
        file("${projectDir}/assets/orthodb_sources.json", checkIfExists: true))

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
    ACQUIRE_PROTEIN_SOURCE(routed.acquire, proteinHelper)
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
    COMBINE_PROTEIN_EVIDENCE(combinations, proteinHelper)
    REPORT_PROTEIN_EVIDENCE(COMBINE_PROTEIN_EVIDENCE.out.evidence.map { aid, bundle -> bundle }.collect(),
        file("${projectDir}/py_scripts/report_protein_evidence.py", checkIfExists: true))
    ASSERT_PROTEIN_EVIDENCE(REPORT_PROTEIN_EVIDENCE.out.status)
    // Assembly-independent RNA preparation and raw/filtered QC.
    PLAN_RNA_PREPARATION(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out,
        encode([quality: params.rna_trim_quality, window: params.rna_trim_window,
                min_length: params.rna_min_read_length, max_n: params.rna_max_n_bases]), rnaHelper)

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
    ACQUIRE_RNA_READSET(acquisition, rnaHelper)
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
        tuple(rid, encode([readset_id: rid, stage: stage]), suffixes, reads[0], reads.size() == 2 ? reads[1] : [])
    }
    QC_RNA_RAW(rawReads.map { rid, manifest, reads -> qcInput('raw', rid, manifest, reads) }, rnaHelper)
    def trimming = rawReads.join(readsets).map { rid, manifest, reads, row ->
        tuple(rid, encode([readset_id: rid] + row.trim_request), manifest, reads)
    }
    TRIM_RNA_READSET(trimming, rnaHelper)
    def clean = splitStatus(TRIM_RNA_READSET.out.bundle)
    def resolveReads = { rid, bundle, priorManifest, priorReads ->
        def metadata = new groovy.json.JsonSlurper().parse(new File("${bundle}/rna.json"))
        tuple(rid, file("${bundle}/rna.json"), metadata.passthrough ? priorReads :
            (metadata.status == 'ready' ? metadata.files.collect { file("${bundle}/${it}") } : []))
    }
    def cleanReads = clean.ready.join(rawReads).map(resolveReads)
    QC_RNA_PROCESSED(cleanReads.map { rid, manifest, reads -> qcInput('processed', rid, manifest, reads) }, rnaHelper)
    def selecting = cleanReads.join(readsets).map { rid, manifest, reads, row ->
        tuple(rid, encode([readset_id: rid] + row.subset_request), manifest, reads)
    }
    SUBSET_RNA_READSET(selecting, rnaHelper)
    def selectedReads = SUBSET_RNA_READSET.out.bundle.join(cleanReads).map(resolveReads)
    def failedReads = raw.failed.mix(clean.failed).map { rid, bundle -> tuple(rid, file("${bundle}/rna.json"), []) }
    def finalStates = selectedReads.mix(failedReads).toList()
    def finalResults = finalStates.map { rows -> rows.collect { it[1] } }
    def readPaths = finalStates.map { rows -> encode(rows.collectEntries { rid, manifest, reads ->
        [rid, reads.collect { it.toString() }]
    }) }
    // Stage individual ZIPs so MultiQC does not need to traverse linked directories.
    def collectQc = { channel -> channel.flatMap { rid, bundle ->
        def reports = file("${bundle}/*_fastqc.zip", checkIfExists: true)
        reports instanceof Collection ? reports : [reports]
    }.toList() }
    def rawQcResults = collectQc(QC_RNA_RAW.out.qc)
    def processedQcResults = collectQc(QC_RNA_PROCESSED.out.qc)
    def fastpReports = clean.ready.flatMap { rid, bundle ->
        def report = file("${bundle}/${rid}_fastp.json")
        report.exists() ? [report] : []
    }.toList()
    REPORT_RNA_PREPARATION(PLAN_RNA_PREPARATION.out.bundle, finalResults, readPaths, rawQcResults, processedQcResults, fastpReports,
        file("${projectDir}/py_scripts/report_rna_preparation.py"))
    ASSERT_RNA_PREPARATION(REPORT_RNA_PREPARATION.out.status)
    // Reference normalization and existing repeat discovery/curation/masking.
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
    }, repeatHelper)
    def referenceStates = FILTER_ASSEMBLY.out.bundle.map { aid, folder ->
        tuple(aid, file("${folder}/${aid}.reference.json"))
    }
    def referenceSplit = referenceStates.branch {
        ready: parse(it[1]).status == 'ready'
        failed: true
    }
    def preparedReferences = referenceSplit.ready.join(FILTER_ASSEMBLY.out.bundle).map { aid, state, folder ->
        tuple(aid, file("${folder}/${aid}.filtered.fasta"), state)
    }
    def librarySources
    def teLibrarySources
    if (params.repeat_library) {
        // An explicit complete library is a scientific input, also used by tiny masking fixtures.
        def supplied = file(params.repeat_library, checkIfExists: true)
        librarySources = preparedReferences.map { aid, reference, state -> tuple(aid, supplied) }
        teLibrarySources = librarySources
    } else {
        if (!params.dfam_database) error 'Provide --dfam_database for the existing Dfam FamDB directory'
        if (!params.diamond_database && !params.repeat_filter_database) error 'Provide --diamond_database or --repeat_filter_database'
        def dfam = file("${params.dfam_database}/dfam_db", checkIfExists: true)
        def generalDb = file(params.repeat_filter_database ?: "${params.diamond_database}/swissprot.dmnd", checkIfExists: true)
        def requestedTaxa = preparedReferences.join(assemblies).map { aid, reference, state, row ->
            tuple(taxonKey(row.settings.taxa_name), aid, row.settings.taxa_name)
        }
        CURATE_DFAM(requestedTaxa.map { key, aid, taxa -> tuple(key, taxa, dfam) }.unique { it[0] })
        def dfamByAssembly = requestedTaxa.combine(CURATE_DFAM.out.bundle, by: 0).map { key, aid, taxa, folder ->
            tuple(aid, file("${folder}/dfam_subset.lib"), file("${folder}/dfam_subset.source.json"))
        }
        REPEAT_MODEL(preparedReferences.join(dfamByAssembly).map { aid, reference, state, library, source ->
            tuple(aid, reference, source)
        }, params.repeat_seed)
        CURATE_REPEATS(REPEAT_MODEL.out.families_fa, repeatHelper)
        CLASSIFY_REPEATS(CURATE_REPEATS.out.clustered_repeats)
        def classified
        if (params.repeat_reclassify) {
            RECLASSIFY_REPEATS(CLASSIFY_REPEATS.out.curated_repeats, file("${projectDir}/py_scripts/run_terrier.py"))
            classified = RECLASSIFY_REPEATS.out.reclassified_lib
        } else {
            classified = CLASSIFY_REPEATS.out.curated_repeats
        }
        EXTRACT_REPEAT_PEPS()
        teLibrarySources = classified
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
    AUDIT_REPEAT_LIBRARY(librarySources, repeatHelper)
    def acceptedLibraries = AUDIT_REPEAT_LIBRARY.out.status.filter { aid, status -> parse(status).status == 'ready' }
        .join(librarySources).map { aid, status, library -> tuple(aid, library) }
    def maskInputs = preparedReferences.join(acceptedLibraries).map { aid, reference, state, library ->
        tuple(aid, reference, library)
    }
    REPEAT_MASK_PARALLEL(maskInputs, params.n_rm_chunks)
    def benchmarkChecks = Channel.value([])
    if (params.repeat_mask_benchmark) {
        UNSPLIT_MASK_BENCHMARK(maskInputs.map { aid, reference, library -> tuple(aid + '__benchmark', reference, library) }, 1)
        def benchmark = UNSPLIT_MASK_BENCHMARK.out.masked_fasta.join(UNSPLIT_MASK_BENCHMARK.out.repeat_out)
            .map { aid, masked, annotations -> tuple(aid.toString().replaceFirst(/__benchmark$/, ''), masked, annotations) }
        COMPARE_REPEAT_MASKS(REPEAT_MASK_PARALLEL.out.masked_fasta.join(REPEAT_MASK_PARALLEL.out.repeat_out)
            .join(benchmark), repeatHelper)
        benchmarkChecks = COMPARE_REPEAT_MASKS.out.toList()
    }
    def validations = REPEAT_MASK_PARALLEL.out.masked_fasta.join(REPEAT_MASK_PARALLEL.out.repeat_out)
        .join(referenceSplit.ready).join(acceptedLibraries).map { aid, masked, annotations, state, library ->
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
    ]), repeatHelper)
    REPORT_REFERENCE_REPEATS(RESOLVE_ANNOTATION_INPUTS.out.bundle, referenceStates.map { it[1] }.toList(),
        AUDIT_REPEAT_LIBRARY.out.status.map { it[1] }.toList(), VALIDATE_REFERENCE_REPEATS.out.status.map { it[1] }.toList(),
        benchmarkChecks,
        file("${projectDir}/py_scripts/report_reference_repeats.py"))
    ASSERT_REFERENCE_REPEATS(REPORT_REFERENCE_REPEATS.out.status, benchmarkChecks)
    // RNA alignment can run once normalization and RNA preparation are ready, alongside repeats.
    PLAN_CONNECTED_RNA_ALIGNMENT(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out,
        REPORT_RNA_PREPARATION.out.status, ASSERT_RNA_PREPARATION.out,
        referenceSplit.ready.join(FILTER_ASSEMBLY.out.bundle).map { aid, state, folder -> folder }.toList(),
        encode([reference_dir: 'reference_inputs', strand_policy: [min_informative: params.rna_strand_min_informative,
            dominant_fraction: params.rna_strand_dominant_fraction, unstranded_tolerance: params.rna_strand_unstranded_tolerance]]), alignmentHelper)
    def alignmentPlans = PLAN_CONNECTED_RNA_ALIGNMENT.out.plan.map { path -> new groovy.json.JsonSlurper().parse(path.toFile()) }
    def alignmentReferences = alignmentPlans.flatMap { plan -> plan.references }
    INDEX_RNA_REFERENCE(alignmentReferences.map { row ->
        tuple(row.reference_id, file(row.reference, checkIfExists: true), file(row.identity, checkIfExists: true))
    }.unique { it[0] }, helpers)
    def alignmentJobs = alignmentPlans.flatMap { plan -> plan.jobs }.map { job -> tuple(job.reference_id, job) }
    def indexed = alignmentJobs.combine(INDEX_RNA_REFERENCE.out.index, by: 0).map { refid, job, index ->
        tuple([job.assembly_id, job.readset_id], job, index)
    }
    ALIGN_RNA_READSET(indexed.map { key, job, index ->
        tuple(job.assembly_id, job.readset_id, encode(job.spec), index,
              job.read_paths.collect { file(it, checkIfExists: true) })
    }, readHelpers)
    def aligned = ALIGN_RNA_READSET.out.alignment.map { aid, rid, folder -> tuple([aid, rid], folder) }
    QC_RNA_ALIGNMENT(indexed.join(aligned).map { key, job, index, folder ->
        tuple(job.assembly_id, job.readset_id, encode(job.spec + [min_mapq: params.rna_strand_min_mapq]), folder, index)
    }, alignmentHelper)
    def qc = QC_RNA_ALIGNMENT.out.qc.map { aid, rid, folder -> tuple([aid, rid], folder) }
    RESOLVE_RNA_STRAND(indexed.join(qc).map { key, job, index, folder ->
        tuple(job.assembly_id, job.readset_id, encode(job.policy), folder)
    }, alignmentHelper)
    def alignmentBamLocations = ALIGN_RNA_READSET.out.alignment.map { aid, rid, folder ->
        [assembly_id: aid, readset_id: rid, bam: "readsets/${aid}/${rid}/aligned.bam",
         bam_index: "readsets/${aid}/${rid}/aligned.bam.csi"]
    }.toList().map { rows -> encode(rows) }
    def alignmentLogs = ALIGN_RNA_READSET.out.alignment.map { aid, rid, folder ->
        [aid: aid, rid: rid, folder: folder]
    }.toList().map { rows ->
        def samples = []
        def logs = []
        rows.sort { a, b -> (a.aid <=> b.aid) ?: (a.rid <=> b.rid) }.each { row ->
            ['hisat2.log', 'flagstat.txt', 'samtools_stats.txt'].each { name ->
                samples << "${row.aid}/${row.rid}".toString()
                logs << file("${row.folder}/${name}", checkIfExists: true)
            }
        }
        tuple(samples, logs)
    }
    REPORT_RNA_ALIGNMENT(PLAN_CONNECTED_RNA_ALIGNMENT.out.plan,
        QC_RNA_ALIGNMENT.out.qc.map { aid, rid, folder -> file("${folder}/alignment.json", checkIfExists: true) }.toList(),
        RESOLVE_RNA_STRAND.out.strand.map { aid, rid, path -> path }.toList(), alignmentBamLocations, alignmentLogs,
        [alignmentHelper, file("${projectDir}/py_scripts/report_rna_alignment.py", checkIfExists: true)])
    ASSERT_RNA_ALIGNMENT(REPORT_RNA_ALIGNMENT.out.status)
    // Independent genome tracks reuse native reference channels; no prediction dependency.
    if (!params.braker_check_only) GENOME_TRACKS(preparedReferences.join(assemblies).join(FILTER_ASSEMBLY.out.bundle).map { aid,genome,state,row,folder ->
        tuple(aid,row.settings,genome,file("${folder}/${aid}.contigs.json",checkIfExists:true))
    },REPEAT_MASK_PARALLEL.out.masked_fasta)
    // Gene prediction: join scientific files by assembly, never read asynchronous publishDir outputs.
    PLAN_BRAKER_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle, REPORT_REFERENCE_REPEATS.out.status,
        REPORT_PROTEIN_EVIDENCE.out.status, REPORT_RNA_ALIGNMENT.out.status, brakerHelpers)
    def predictionRequests = PLAN_BRAKER_INPUTS.out.plan.flatMap { path -> parse(path).jobs }
        .map { row -> tuple(row.assembly_id, encode(row), row) }
    // Gather BAM tuples as optional evidence: an inner join would drop RNA-free assemblies.
    def nativeAlignments = ALIGN_RNA_READSET.out.alignment.toList().map { rows -> [rows] }
    def predictionFiles = REPEAT_MASK_PARALLEL.out.masked_fasta
        .join(COMBINE_PROTEIN_EVIDENCE.out.evidence).join(FILTER_ASSEMBLY.out.bundle)
        .combine(nativeAlignments)
    def predictionInputs = predictionRequests.join(predictionFiles).map { aid, spec, row, genome, proteins, reference, alignmentRows ->
        def selected = row.readsets.collect { readset ->
            def found = alignmentRows.find { item -> item[0] == aid && item[1] == readset.readset_id }
            if (!found) error "Missing native BAM for ${aid}/${readset.readset_id}"
            found[2]
        }
        tuple(aid, spec, genome,
            row.mode in ['protein', 'rna_protein'] ? file("${proteins}/proteins.fa", checkIfExists: true) : [],
            selected.collect { file("${it}/aligned.bam", checkIfExists: true) },
            selected.collect { file("${it}/aligned.bam.csi", checkIfExists: true) },
            file("${reference}/${aid}.contigs.json", checkIfExists: true))
    }
    VALIDATE_BRAKER_INPUTS(predictionInputs, brakerHelpers)
    CHECK_BRAKER_SOFTWARE(predictionInputs.join(VALIDATE_BRAKER_INPUTS.out.status), brakerHelpers)
    def predictionResults = Channel.value([])
    if (!params.braker_check_only) {
        BRAKER(predictionInputs.join(CHECK_BRAKER_SOFTWARE.out.runtime), brakerHelpers)
        VALIDATE_BRAKER_OUTPUTS(predictionInputs.join(BRAKER.out.prediction).map { aid, spec, genome, proteins, bams, indices, identity, prediction ->
            tuple(aid, spec, prediction, identity)
        }, brakerHelpers)
        predictionResults = VALIDATE_BRAKER_OUTPUTS.out.status.map { aid, path -> path }.toList()
        def finalSettings = assemblies.map { aid,row -> tuple(aid,
            [run_busco: row.settings.run_busco, run_omark: row.settings.run_omark,
             busco_lineage: row.settings.busco_lineage, taxid: row.settings.taxid]) }
        def postprocessing = predictionInputs.join(BRAKER.out.prediction).join(VALIDATE_BRAKER_OUTPUTS.out.status)
            .join(finalSettings).map { aid,spec,genome,proteins,bams,indices,identity,prediction,validation,settings ->
                def metadata = new groovy.json.JsonSlurperClassic().parseText(new String(spec.decodeBase64(),'UTF-8'))
                tuple(aid,settings,metadata.mode,metadata.readsets,genome,prediction,validation,identity,bams,indices)
            }
        def proteinAlignmentInputs = REPEAT_MASK_PARALLEL.out.masked_fasta.join(COMBINE_PROTEIN_EVIDENCE.out.evidence)
            .join(FILTER_ASSEMBLY.out.bundle).join(assemblies).filter { aid,genome,proteins,reference,row ->
                row.prediction_mode in ['protein','rna_protein']
            }.map { aid,genome,proteins,reference,row ->
                tuple(aid,genome,file("${proteins}/proteins.fa"),file("${reference}/${aid}.contigs.json"),file("${proteins}/evidence.json"))
            }
        ANNOTATION_POSTPROCESSING(postprocessing,proteinAlignmentInputs,
            REPEAT_MASK_PARALLEL.out.repeat_gff.join(teLibrarySources))
        if (params.run_functional_annotation) {
            FUNCTIONAL_ANNOTATION(assemblies.join(ANNOTATION_POSTPROCESSING.out.identified_annotation).map { aid,row,folder ->
                tuple(aid,row.settings,folder)
            },ANNOTATION_POSTPROCESSING.out.omark_classes)
        }
    }
    REPORT_BRAKER(PLAN_BRAKER_INPUTS.out.plan, predictionResults,
        CHECK_BRAKER_SOFTWARE.out.runtime.map { aid, folder -> file("${folder}/inputs.json") }.toList(),
        params.braker_check_only, '../rna_alignment/multiqc/multiqc_report.html', file("${projectDir}/py_scripts/report_braker.py"))
    ASSERT_BRAKER(REPORT_BRAKER.out.status)
}
