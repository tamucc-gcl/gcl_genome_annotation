// Native scientific modules shared with the accepted development entries.
include { PREPARE_UTR_INPUTS; RUN_UTR_REFINEMENT; FINALIZE_UTR_ANNOTATION; EXTRACT_UTR_SEQUENCES; VALIDATE_UTR_REFINEMENT } from '../modules/utrRefinement.nf'
include { ASSEMBLE_AED_RNA; MERGE_AED_RNA; ALIGN_AED_PROTEINS } from '../modules/evidenceAed.nf'
include { FINALIZE_ANNOTATION; REPORT_FINAL_AED } from '../modules/finalAnnotation.nf'
include { PREPARE_INGENANNOT_EVIDENCE; RUN_INGENANNOT_AED } from '../modules/ingenannotAed.nf'
include { ANNOTATION_BASELINE as FINAL_ANNOTATION_QC } from './annotationBaseline.nf'
include { JUNCTION_SUPPORT } from './junctionSupport.nf'
include { MODEL_ASSESSMENT } from './modelAssessment.nf'
include { MODEL_SELECTION } from './modelSelection.nf'
include { FINAL_STRUCTURAL_ANNOTATION } from './finalStructuralAnnotation.nf'

workflow ANNOTATION_POSTPROCESSING {
    take:
    // aid, QC settings, mode, readsets, genome, BRAKER directory, validation,
    // reference identity, native BAMs, native BAM indices
    jobs
    proteinRequests // Independent of prediction: alignment can overlap BRAKER.
    repeatRequests // Native RepeatMasker GFF and pre-filter repeat library for assessment.

    main:
    if (params.run_finalization && !params.run_model_selection)
        error 'Final identifiers and structural QC require run_model_selection=true'
    if (params.run_model_selection && !params.run_model_assessment)
        error 'Model selection requires run_model_assessment=true'
    if (params.ingenannot_rna_penalty < 0 || params.ingenannot_rna_penalty > 1)
        error 'ingenannot_rna_penalty must be between 0 and 1'
    utrHelpers = ['utr_refinement.py','annotation_models.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
    }
    finalHelpers = ['final_annotation.py','ingenannot_aed.py','ingenannot_inputs.py','annotation_models.py','prepare_protein_evidence.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
    }
    rnaJobs = jobs.map { aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices ->
        if (!aid || !prediction || !(mode in ['rna_protein','rna','protein','ab_initio']))
            error "Invalid postprocessing job: assembly=${aid}, mode=${mode}, prediction=${prediction}"
        if (mode in ['rna','rna_protein'] && !readsets) error "Missing RNA readsets for ${aid}"
        def retained = file("${prediction}/braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff")
        tuple(aid, mode in ['rna','rna_protein'] ? readsets : [], bams, indices,
            mode in ['rna','rna_protein'] && !params.aed_reassemble_rna && retained.exists() ? [retained] : [])
    }
    ASSEMBLE_AED_RNA(rnaJobs.filter { it[1] && !it[4] }.flatMap { aid,readsets,bams,indices,retained ->
        readsets.withIndex().collect { row,i -> tuple(aid,row.readset_id,row.strandedness,bams[i],indices[i]) }
    })
    MERGE_AED_RNA(ASSEMBLE_AED_RNA.out.transcripts.groupTuple())
    rnaSources = rnaJobs.filter { !it[1] || it[4] }.map { aid,readsets,bams,indices,retained -> tuple(aid,retained) }
        .mix(MERGE_AED_RNA.out.transcripts.map { aid,path -> tuple(aid,[path]) })
    utrRequests = jobs.join(rnaSources).map { aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,rna ->
        def enabled = params.utr_refine && mode in ['rna','rna_protein']
        tuple(aid,enabled,enabled ? readsets : [],file("${prediction}/${aid}.gff3", checkIfExists: true),
            validation,identity,enabled ? rna : [],enabled ? bams : [],enabled ? indices : [])
    }
    PREPARE_UTR_INPUTS(utrRequests,utrHelpers)
    RUN_UTR_REFINEMENT(utrRequests.join(PREPARE_UTR_INPUTS.out.status.map { aid,path ->
        // Invoke the parser directly: workflow-bound closures are not script
        // methods inside operator callbacks on Nextflow 23.10.1.
        tuple(aid,new groovy.json.JsonSlurperClassic().parse(new File(path.toString())))
    }))
    FINALIZE_UTR_ANNOTATION(utrRequests.map { aid,enabled,readsets,gff,validation,identity,rna,bams,indices ->
        tuple(aid,!enabled,gff,identity)
    }.join(RUN_UTR_REFINEMENT.out.annotation).map { aid,skipped,original,identity,raw ->
        tuple(aid,skipped,original,raw,identity)
    },utrHelpers)
    genomes = jobs.map { aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices -> tuple(aid,genome) }
    modelFiles = utrRequests.map { aid,enabled,readsets,gff,validation,identity,rna,bams,indices -> tuple(aid,!enabled,gff,identity) }
        .join(genomes).join(FINALIZE_UTR_ANNOTATION.out.annotation)
    EXTRACT_UTR_SEQUENCES(modelFiles.map { aid,skipped,original,identity,genome,refined -> tuple(aid,genome,original,refined) })
    VALIDATE_UTR_REFINEMENT(modelFiles.join(PREPARE_UTR_INPUTS.out.status).join(EXTRACT_UTR_SEQUENCES.out.sequences)
        .join(FINALIZE_UTR_ANNOTATION.out.rejections).map {
            aid,skipped,original,identity,genome,refined,inputs,sequences,rejections ->
            tuple(aid,skipped,original,refined,identity,genome,inputs,sequences,rejections)
        },utrHelpers)
    FINALIZE_ANNOTATION(jobs.join(VALIDATE_UTR_REFINEMENT.out.annotation).join(VALIDATE_UTR_REFINEMENT.out.status).map {
        aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,gff,utrStatus ->
        tuple(aid,mode,['gff3','proteins.fa','cds.fa'].collect { suffix ->
            file(suffix == 'gff3' ? "${prediction}/${aid}.gff3" : "${prediction}/${aid}_${suffix}", checkIfExists: true)
        },validation,gff,utrStatus)
    },['final_annotation_inputs.py','prepare_protein_evidence.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
    })
    junctionReports = Channel.empty()
    aedReports = Channel.empty()
    selectedAnnotations = Channel.empty()
    identifiedAnnotations = Channel.empty()
    structuralReports = Channel.empty()
    selectedOmarkClasses = Channel.empty()
    beforeQc = Channel.empty()
    if (params.run_junction_qc) {
        JUNCTION_SUPPORT(jobs.map { aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices ->
            tuple(aid,readsets,genome,identity,bams,indices)
        },FINALIZE_ANNOTATION.out.bundle.map { aid,folder ->
            tuple(aid,file("${folder}/${aid}.gff3"),file("${folder}/annotation_status.json"))
        })
        junctionReports = JUNCTION_SUPPORT.out.transcripts.join(JUNCTION_SUPPORT.out.status)
    }
    if (params.run_baseline_qc) {
        FINAL_ANNOTATION_QC(jobs.join(FINALIZE_ANNOTATION.out.bundle).map {
            aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,folder ->
            tuple(aid,settings,genome,file("${folder}/${aid}.gff3"),file("${folder}/${aid}_proteins.fa"),file("${folder}/annotation_status.json"))
        },['annotation_baseline.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
            file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
        })
        beforeQc = FINAL_ANNOTATION_QC.out.agat.mix(FINAL_ANNOTATION_QC.out.busco).mix(FINAL_ANNOTATION_QC.out.omark)
            .mix(FINAL_ANNOTATION_QC.out.selections.map { aid,folder -> tuple(aid,file("${folder}/selection.json")) })
    }
    if (params.run_aed) {
        ALIGN_AED_PROTEINS(proteinRequests,['prepare_protein_evidence.py','reference_repeats.py'].collect {
            file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
        })
        proteinSources = ALIGN_AED_PROTEINS.out.alignment.toList().map { rows -> [rows] }
        aedRequests = jobs.join(FINALIZE_ANNOTATION.out.bundle).join(rnaSources).combine(proteinSources).map {
            aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,folder,rna,proteinRows ->
            def match = proteinRows.find { it[0] == aid }
            if (mode in ['protein','rna_protein'] && !match) error "Missing protein alignment for ${aid}"
            tuple(aid,file("${folder}/${aid}.gff3"),file("${folder}/annotation_status.json"),rna,match ? [match[1]] : [])
        }
        PREPARE_INGENANNOT_EVIDENCE(aedRequests,['ingenannot_inputs.py','prepare_protein_evidence.py'].collect {
            file("${params.pipeline_root}/py_scripts/${it}", checkIfExists: true)
        })
        aedReady = aedRequests.map { aid,gff,validation,rna,protein -> tuple(aid,gff) }
            .join(PREPARE_INGENANNOT_EVIDENCE.out.evidence)
        RUN_INGENANNOT_AED(aedReady)
        REPORT_FINAL_AED(aedReady.join(RUN_INGENANNOT_AED.out.annotated).join(jobs.map {
            aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices -> tuple(aid,identity)
        }).map { aid,gff,indexed,annotated,identity -> tuple(aid,gff,identity,indexed,annotated) },finalHelpers)
        aedReports = REPORT_FINAL_AED.out.scores.join(REPORT_FINAL_AED.out.status)
    }
    if (params.run_model_assessment) {
        MODEL_ASSESSMENT(jobs.join(FINALIZE_ANNOTATION.out.bundle).map {
            aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,folder ->
            // Coding sequences are unchanged by validated UTR refinement.
            // Keep their original native paths so UTR/report changes do not
            // invalidate the TE protein search through newly published copies.
            tuple(aid,genome,file("${folder}/${aid}.gff3"),file("${prediction}/${aid}_proteins.fa"),
                file("${prediction}/${aid}_cds.fa"),file("${folder}/annotation_status.json"),identity)
        },repeatRequests,aedReports,junctionReports)
        if (params.run_model_selection) {
            MODEL_SELECTION(jobs.join(FINALIZE_ANNOTATION.out.bundle).map {
                aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices,folder ->
                tuple(aid,file("${folder}/${aid}.gff3"),file("${folder}/${aid}_proteins.fa"),
                    file("${folder}/${aid}_cds.fa"),file("${folder}/annotation_status.json"),identity)
            },MODEL_ASSESSMENT.out.transcripts.join(MODEL_ASSESSMENT.out.status))
            selectedAnnotations = MODEL_SELECTION.out.annotation
            if (params.run_finalization) {
                FINAL_STRUCTURAL_ANNOTATION(jobs.map { aid,settings,mode,readsets,genome,prediction,validation,identity,bams,indices ->
                    tuple(aid,settings,genome,identity)
                },selectedAnnotations,MODEL_ASSESSMENT.out.transcripts.join(MODEL_ASSESSMENT.out.status),beforeQc)
                identifiedAnnotations = FINAL_STRUCTURAL_ANNOTATION.out.annotation
                structuralReports = FINAL_STRUCTURAL_ANNOTATION.out.report
                selectedOmarkClasses = FINAL_STRUCTURAL_ANNOTATION.out.omark_classes
            }
        }
    }

    emit:
    annotation = FINALIZE_ANNOTATION.out.bundle
    selected_annotation = selectedAnnotations
    identified_annotation = identifiedAnnotations
    omark_classes = selectedOmarkClasses
    utr_status = VALIDATE_UTR_REFINEMENT.out.status
    structural_reports = structuralReports
    selection_status = selectedAnnotations.map { aid,folder -> tuple(aid,file("${folder}/selection_status.json",checkIfExists:true)) }
}
