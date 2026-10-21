include { RUN_PORTCULLIS; REPORT_JUNCTION_SUPPORT } from '../modules/junctionSupport.nf'

workflow JUNCTION_SUPPORT {
    take:
    // aid, readsets, genome, identity, BAMs, CSI indices. No annotation input:
    // annotation/report edits must not invalidate native junction detection.
    alignments
    annotations // aid, final GFF, validated final status

    main:
    if (params.portcullis_threshold < 0 || params.portcullis_threshold > 1)
        error 'portcullis_threshold must be between 0 and 1'
    if (!(params.portcullis_paired_orientation in ['FR','RF','FF','UNKNOWN']))
        error 'portcullis_paired_orientation must be FR, RF, FF or UNKNOWN'
    RUN_PORTCULLIS(alignments.flatMap { aid,readsets,genome,identity,bams,indices ->
        if (readsets.size() != bams.size() || bams.size() != indices.size())
            error "Portcullis readset/BAM count mismatch for ${aid}"
        if (readsets.collect { it.readset_id }.toSet().size() != readsets.size())
            error "Duplicate Portcullis readset for ${aid}"
        if (!(aid ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || readsets.any { !(it.readset_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) })
            error 'Portcullis assembly/readset IDs must be safe file identifiers'
        readsets.withIndex().collect { row,i ->
            def strand = [reverse:'firststrand',forward:'secondstrand',unstranded:'unstranded'][row.strandedness]
            if (!strand) error "Unresolved Portcullis library strand: ${aid}/${row.readset_id}"
            def orientation = row.layout == 'single' ? 'SE' : params.portcullis_paired_orientation
            def request = groovy.json.JsonOutput.toJson([assembly_id:aid,readset:row,orientation:orientation]).getBytes('UTF-8').encodeBase64().toString()
            tuple(aid,row.readset_id,request,strand,orientation,genome,identity,bams[i],indices[i])
        }
    },['portcullis_inputs.py','reference_repeats.py','prepare_protein_evidence.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists:true)
    })
    // Optional evidence must not drop RNA-free assemblies in an inner join.
    nativeRows = RUN_PORTCULLIS.out.junctions.toList().map { rows -> [rows] }
    REPORT_JUNCTION_SUPPORT(annotations.join(alignments.map {
        aid,readsets,genome,identity,bams,indices -> tuple(aid,readsets,identity)
    }).combine(nativeRows).map { aid,gff,validation,readsets,identity,rows ->
        def selected = rows.findAll { it[0] == aid }.sort { a,b -> a[1] <=> b[1] }
        def request = groovy.json.JsonOutput.toJson([assembly_id:aid,readsets:readsets]).getBytes('UTF-8').encodeBase64().toString()
        tuple(aid,request,gff,validation,identity,selected.collect { it[2] })
    },['junction_support.py','portcullis_inputs.py','annotation_models.py','reference_repeats.py','prepare_protein_evidence.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}", checkIfExists:true)
    })
    emit:
    status = REPORT_JUNCTION_SUPPORT.out.status
    transcripts = REPORT_JUNCTION_SUPPORT.out.transcripts
}
