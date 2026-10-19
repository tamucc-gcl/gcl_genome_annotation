// Copied to the private pipeline root by validate_postprocessing_cluster.py.
nextflow.enable.dsl = 2
include { ANNOTATION_POSTPROCESSING } from './workflows/annotationPostprocessing.nf'

workflow {
    def root = file(params.utr_inputs_root).toString()
    // The two channels read the same rows concurrently. Use eager ordinary maps;
    // JsonSlurper's LazyMap can expose an incompletely built map on another thread.
    def rows = new groovy.json.JsonSlurperClassic().parse(new File("${root}/fixtures.json"))
    def jobs = Channel.fromList(rows).map { row ->
        def aid = row.assembly_id
        if (!aid || !(row.mode in ['rna_protein','rna','protein','ab_initio']))
            error "Invalid postprocessing fixture metadata: ${row}"
        tuple(aid,[run_busco:false,run_omark:false,busco_lineage:null,taxid:null],row.mode,row.readsets,
            file("${root}/repeats/${aid}.masked.fasta"),file("${root}/structural_annotation/braker/${aid}"),
            file("${root}/structural_annotation/braker_status/${aid}.prediction.json"),
            file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json"),
            row.readsets.collect { file("${root}/rna_alignment/${it.bam}") },
            row.readsets.collect { file("${root}/rna_alignment/${it.bam_index}") })
    }
    def proteins = Channel.fromList(rows).filter { it.mode in ['protein','rna_protein'] }.map { row ->
        def aid = row.assembly_id
        if (!aid) error "Missing assembly_id in postprocessing protein fixture: ${row}"
        tuple(aid,file("${root}/repeats/${aid}.masked.fasta"),file("${root}/protein_evidence/${aid}/proteins.fa"),
            file("${root}/filtered_assemblies/${aid}/${aid}.contigs.json"),file("${root}/protein_evidence/${aid}/evidence.json"))
    }
    ANNOTATION_POSTPROCESSING(jobs,proteins,Channel.empty())
}
