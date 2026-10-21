#!/usr/bin/env nextflow
// Development/pilot entry. Normal main.nf passes native task outputs directly.
nextflow.enable.dsl = 2
include { FUNCTIONAL_ANNOTATION } from '../../workflows/functionalAnnotation.nf'

workflow {
    if (!params.functional_inputs_root) error 'Provide --functional_inputs_root with accepted final identified outputs'
    def root = file(params.functional_inputs_root).toAbsolutePath().toString()
    def state = new groovy.json.JsonSlurperClassic().parse(new File("${root}/preflight/resolved_inputs.json"))
    if (state.status != 'ready') error 'Accepted preflight is not ready'
    def rows = state.assemblies
    def jobs = Channel.fromList(rows).map { row ->
        def settings = [:]+row.settings
        // Global switches may disable an accepted per-assembly analysis for a pilot.
        ['interproscan','eggnog','diamond','kofamscan','deepsig'].each { tool ->
            settings['run_'+tool] = settings['run_'+tool] && params['run_'+tool]
        }
        if (params.eggnog_taxid != null) settings.eggnog_taxid = params.eggnog_taxid.toString()
        tuple(row.assembly_id,settings,file("${root}/structural_annotation/identified/${row.assembly_id}",checkIfExists:true))
    }
    def omark = Channel.fromList(rows).filter { it.settings.run_omark }.map { row ->
        tuple(row.assembly_id,file("${root}/qc/selected/${row.assembly_id}/omark/${row.assembly_id}.ump",checkIfExists:true))
    }
    if (params.run_functional_annotation) FUNCTIONAL_ANNOTATION(jobs,omark)
}
