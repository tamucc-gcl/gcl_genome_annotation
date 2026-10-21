include { PACKAGE_ANNOTATION_RELEASE; ANNOTATION_RELEASE_INDEX } from '../modules/annotationRelease.nf'

workflow ANNOTATION_RELEASE {
    take:
    references // aid, unmasked genome, masked genome, identity, repeats, repeat validation
    identified // aid, accepted identified bundle
    functional // aid, decorated coding GFF, table, gene products, functional validation; empty if disabled
    tracks // aid, track, selected GFF
    trackStatuses // aid, complete track manifest (including disabled tracks)
    main:
    optionalRows=functional.toList().map { [it] }
    trackRows=tracks.toList().map { [it] }
    PACKAGE_ANNOTATION_RELEASE(references.join(identified).join(trackStatuses).combine(optionalRows).combine(trackRows).map {
        aid,genome,masked,identity,repeats,repeatStatus,folder,trackStatus,functions,selected ->
        def matches=functions.findAll { it[0]==aid }
        if (matches.size()!=(params.run_functional_annotation ? 1 : 0)) error "Missing/duplicate/unexpected functional release input for ${aid}"
        def f=matches ? matches[0] : [aid,[],[],[],[]]
        tuple(aid,folder,genome,masked,identity,repeats,repeatStatus,trackStatus,
            selected.findAll { it[0]==aid }.sort { it[1] }.collect { it[2] },params.run_functional_annotation,f[1],f[2],f[3],f[4])
    },['package_annotation.py','functional_report.py','functional_inputs.py','model_metrics.py','annotation_models.py',
        'prepare_protein_evidence.py','reference_repeats.py'].collect { file("${params.pipeline_root}/py_scripts/${it}",checkIfExists:true) })
    releases=PACKAGE_ANNOTATION_RELEASE.out.bundle.toList().map { it.sort { row -> row[0] } }
    // Do not silently advertise a partial release if a required join lost an assembly.
    ANNOTATION_RELEASE_INDEX(references.map { it[0] }.toList().map { it.sort() },releases.map { it.collect { row -> row[1] } })
    emit:
    annotation=PACKAGE_ANNOTATION_RELEASE.out.bundle
}
