#!/usr/bin/env nextflow
// Development pilot: reuse accepted predictions, RNA and BAMs; no BRAKER/repeats.
nextflow.enable.dsl = 2
include { PREPARE_UTR_INPUTS; RUN_UTR_REFINEMENT; FINALIZE_UTR_ANNOTATION; EXTRACT_UTR_SEQUENCES; VALIDATE_UTR_REFINEMENT } from './modules/utrRefinement.nf'
include { ASSEMBLE_AED_RNA; MERGE_AED_RNA } from './modules/evidenceAed.nf'
include { AGAT as UTR_AGAT } from './modules/agat.nf'

workflow {
    def root = file(params.utr_inputs_root ?: params.outdir).toString()
    def parse = { path -> new groovy.json.JsonSlurper().parse(new File(path.toString())) }
    def predictions = parse("${root}/structural_annotation/prediction_status.json")
    if (predictions.status != 'ready') error 'UTR refinement requires accepted BRAKER predictions'
    def helpers = ['utr_refinement.py','annotation_models.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${projectDir}/py_scripts/${it}", checkIfExists: true)
    }
    def base = Channel.fromList(predictions.assemblies).map { row ->
        def aid = row.assembly_id
        if (row.status != 'ready') error "Prediction not ready: ${aid}"
        def enabled = params.utr_refine && row.mode in ['rna','rna_protein']
        def readsets = enabled ? row.readsets : []
        if (enabled && !readsets) error "Missing accepted RNA readsets for ${aid}"
        def candidates = ["${root}/structural_annotation/braker/${aid}/braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff",
            "${root}/structural_annotation/evidence_aed/${aid}/rna_assembly/transcripts_merged.gtf"]
        def retained = enabled ? candidates.find { file(it).exists() } : null
        tuple(aid, enabled, readsets,
            file("${root}/structural_annotation/braker/${aid}/${aid}.gff3", checkIfExists: true),
            file("${root}/structural_annotation/braker_status/${aid}.prediction.json", checkIfExists: true),
            file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json", checkIfExists: true),
            retained ? [file(retained)] : [],
            readsets.collect { file("${root}/rna_alignment/${it.bam}", checkIfExists: true) },
            readsets.collect { file("${root}/rna_alignment/${it.bam_index}", checkIfExists: true) })
    }
    def fallback = base.filter { it[1] && !it[6] }
    ASSEMBLE_AED_RNA(fallback.flatMap { aid, enabled, readsets, gff, validation, identity, rna, bams, indices ->
        readsets.withIndex().collect { row, i -> tuple(aid, row.readset_id, row.strandedness, bams[i], indices[i]) }
    })
    MERGE_AED_RNA(ASSEMBLE_AED_RNA.out.transcripts.groupTuple())
    def rnaSources = base.filter { !it[1] || it[6] }.map { tuple(it[0],it[6]) }.mix(
        MERGE_AED_RNA.out.transcripts.map { aid,path -> tuple(aid,[path]) })
    def requests = base.map { aid, enabled, readsets, gff, validation, identity, rna, bams, indices ->
        tuple(aid,enabled,readsets,gff,validation,identity,bams,indices)
    }.join(rnaSources).map { aid,enabled,readsets,gff,validation,identity,bams,indices,rna ->
        tuple(aid,enabled,readsets,gff,validation,identity,rna,bams,indices)
    }
    PREPARE_UTR_INPUTS(requests,helpers)
    // Content of the approval is a value, not its newly created work path.
    // Reporting/helper edits with identical approved inputs retain ranking cache.
    RUN_UTR_REFINEMENT(requests.join(PREPARE_UTR_INPUTS.out.status.map { aid,path -> tuple(aid,parse(path)) }))
    FINALIZE_UTR_ANNOTATION(requests.map { aid,enabled,readsets,gff,validation,identity,rna,bams,indices ->
        tuple(aid,!enabled,gff,identity)
    }.join(RUN_UTR_REFINEMENT.out.annotation).map { aid,skipped,original,identity,raw_gff ->
        tuple(aid,skipped,original,raw_gff,identity)
    },helpers)
    def genomes = base.map { row -> tuple(row[0],file("${root}/repeats/${row[0]}.masked.fasta", checkIfExists: true)) }
    def modelFiles = requests.map { aid,enabled,readsets,gff,validation,identity,rna,bams,indices -> tuple(aid,!enabled,gff,identity) }
        .join(genomes).join(FINALIZE_UTR_ANNOTATION.out.annotation)
    EXTRACT_UTR_SEQUENCES(modelFiles.map { aid,skipped,original,identity,genome,refined -> tuple(aid,genome,original,refined) })
    VALIDATE_UTR_REFINEMENT(modelFiles.join(PREPARE_UTR_INPUTS.out.status).join(EXTRACT_UTR_SEQUENCES.out.sequences).join(FINALIZE_UTR_ANNOTATION.out.rejections).map {
        aid,skipped,original,identity,genome,refined,inputs,sequences,rejections -> tuple(aid,skipped,original,refined,identity,genome,inputs,sequences,rejections)
    },helpers)
    UTR_AGAT(VALIDATE_UTR_REFINEMENT.out.annotation.join(genomes).map { aid,gff,genome -> tuple(aid,genome,gff) })
}
