#!/usr/bin/env nextflow
// QC-only entry: accepted final models and existing alignments; no prediction.
nextflow.enable.dsl = 2
include { JUNCTION_SUPPORT } from '../../workflows/junctionSupport.nf'

workflow {
    def root = file(params.junction_inputs_root ?: params.outdir).toString()
    def predictions = new groovy.json.JsonSlurperClassic().parse(new File("${root}/structural_annotation/prediction_status.json"))
    if (predictions.status != 'ready') error 'Junction QC requires accepted prediction metadata'
    if (params.run_junction_qc) {
        def rows = predictions.assemblies
        def alignments = Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            if (!aid || row.status != 'ready') error "Invalid junction QC assembly: ${row}"
            def readsets = row.mode in ['rna','rna_protein'] ? row.readsets : []
            if (row.mode in ['rna','rna_protein'] && !readsets) error "Missing RNA evidence: ${aid}"
            tuple(aid,readsets,file("${root}/repeats/${aid}.masked.fasta",checkIfExists:true),
                file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json",checkIfExists:true),
                readsets.collect { file("${root}/rna_alignment/${it.bam}",checkIfExists:true) },
                readsets.collect { file("${root}/rna_alignment/${it.bam_index}",checkIfExists:true) })
        }
        def annotations = Channel.fromList(rows).map { row ->
            def aid = row.assembly_id
            tuple(aid,file("${root}/structural_annotation/final/${aid}/${aid}.gff3",checkIfExists:true),
                file("${root}/structural_annotation/final/${aid}/annotation_status.json",checkIfExists:true))
        }
        JUNCTION_SUPPORT(alignments,annotations)
    }
}
