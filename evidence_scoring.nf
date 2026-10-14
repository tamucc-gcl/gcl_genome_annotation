#!/usr/bin/env nextflow
// Development entry for scoring accepted predictions; no prediction or filtering.
nextflow.enable.dsl = 2
include { ALIGN_AED_PROTEINS; ASSEMBLE_AED_RNA; MERGE_AED_RNA; NORMALIZE_AED_EVIDENCE; SCORE_EVIDENCE_AED } from './modules/evidenceAed.nf'

workflow {
    def root = file(params.aed_inputs_root ?: params.outdir).toString()
    def parse = { path -> new groovy.json.JsonSlurper().parse(new File(path.toString())) }
    def preflight = parse(file("${root}/preflight/resolved_inputs.json", checkIfExists: true))
    def predictions = parse(file("${root}/structural_annotation/prediction_status.json", checkIfExists: true))
    if (preflight.status != 'ready' || predictions.status != 'ready')
        error 'Evidence scoring requires accepted inputs and completed BRAKER predictions'
    def helpers = ['evidence_aed.py', 'prepare_protein_evidence.py'].collect { file("${projectDir}/py_scripts/${it}") }
    def jobs = preflight.assemblies.collect { assembly ->
        def aid = assembly.assembly_id
        def prediction = predictions.assemblies.find { it.assembly_id == aid }
        if (!prediction || prediction.status != 'ready' || prediction.mode != assembly.prediction_mode)
            error "Missing/mismatched prediction for ${aid}"
        [aid: aid, prediction: prediction,
         proteins: assembly.prediction_mode in ['protein', 'rna_protein'],
         rna: assembly.prediction_mode in ['rna', 'rna_protein']]
    }
    def channelJobs = Channel.fromList(jobs)
    ALIGN_AED_PROTEINS(channelJobs.filter { it.proteins }.map { job ->
        tuple(job.aid, file("${root}/repeats/${job.aid}.masked.fasta", checkIfExists: true),
            file("${root}/protein_evidence/assemblies/${job.aid}/proteins.fa", checkIfExists: true),
            file("${root}/filtered_assemblies/${job.aid}/${job.aid}.contigs.json", checkIfExists: true),
            file("${root}/protein_evidence/assemblies/${job.aid}/evidence.json", checkIfExists: true))
    }, ['prepare_protein_evidence.py', 'reference_repeats.py'].collect { file("${projectDir}/py_scripts/${it}") })
    def rnaJobs = channelJobs.filter { it.rna }.map { job ->
        def retained = file("${root}/structural_annotation/braker/${job.aid}/braker_run/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff")
        job + [retained: !params.aed_reassemble_rna && retained.exists() ? retained : null]
    }
    def routed = rnaJobs.branch {
        retained: it.retained != null
        assemble: true
    }
    def existing = routed.retained.map { job -> tuple(job.aid, job.retained) }
    ASSEMBLE_AED_RNA(routed.assemble.flatMap { job ->
        if (!job.prediction.readsets) error "Missing RNA readsets for ${job.aid}"
        job.prediction.readsets.collect { readset ->
            if (!(readset.strandedness in ['reverse', 'forward', 'unstranded']))
                error "Unresolved RNA strand for ${job.aid}/${readset.readset_id}"
            tuple(job.aid, readset.readset_id, readset.strandedness,
                file("${root}/rna_alignment/${readset.bam}", checkIfExists: true),
                file("${root}/rna_alignment/${readset.bam_index}", checkIfExists: true))
        }
    })
    MERGE_AED_RNA(ASSEMBLE_AED_RNA.out.transcripts.groupTuple())
    def rna = existing.mix(MERGE_AED_RNA.out.transcripts).toList()
    def protein = ALIGN_AED_PROTEINS.out.alignment.toList()
    def files = channelJobs.combine(rna.map { [it] }).combine(protein.map { [it] }).map { job, rnaRows, proteinRows ->
        def r = rnaRows.find { it[0] == job.aid }
        def p = proteinRows.find { it[0] == job.aid }
        if (job.rna && !r || job.proteins && !p) error "Missing requested evidence for ${job.aid}"
        tuple(job.aid, file("${root}/filtered_assemblies/${job.aid}/${job.aid}.contigs.json", checkIfExists: true),
            r ? [r[1]] : [], p ? [p[1]] : [])
    }
    NORMALIZE_AED_EVIDENCE(files, helpers)
    SCORE_EVIDENCE_AED(NORMALIZE_AED_EVIDENCE.out.evidence.map { aid, evidence ->
        tuple(aid, file("${root}/structural_annotation/braker/${aid}/${aid}.gff3", checkIfExists: true),
            file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json", checkIfExists: true), evidence,
            file("${root}/structural_annotation/braker_status/${aid}.prediction.json", checkIfExists: true))
    }, helpers + [file("${projectDir}/py_scripts/evidence_aed_selection.py")])
}
