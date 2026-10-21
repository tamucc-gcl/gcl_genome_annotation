include { PREPARE_MODEL_METRICS; REPORT_MODEL_ASSESSMENT } from '../modules/modelAssessment.nf'
include { BUILD_TE_REPEAT_DB; SEARCH_TE_PROTEINS; FLAG_TE_GENES } from '../modules/flagTeGenes.nf'

def validateTeSettings(overlap, evalue) {
    // CLI scientific notation can remain a String in Nextflow 23.10.1.
    def overlapValue
    def evalueValue
    try {
        overlapValue = new BigDecimal(overlap.toString())
        evalueValue = new BigDecimal(evalue.toString())
    } catch (IllegalArgumentException | NullPointerException ex) {
        throw new IllegalArgumentException('TE overlap fraction/e-value must be finite numbers')
    }
    if (overlapValue <= 0 || overlapValue > 1 || evalueValue <= 0)
        throw new IllegalArgumentException('TE overlap fraction must be >0 and <=1; e-value must be >0')
}

workflow MODEL_ASSESSMENT {
    take:
    annotations // aid, genome, final GFF, proteins, CDS, validation, identity
    repeatRequests // aid, RepeatMasker GFF, pre-filter nucleotide repeat library
    aedReports // aid, scores, status
    junctionReports // aid, transcript support, status

    main:
    PREPARE_MODEL_METRICS(annotations,['model_metrics.py','annotation_models.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${params.pipeline_root}/py_scripts/${it}",checkIfExists:true)
    })
    teReports = Channel.value([])
    if (params.run_te_flagging) {
        validateTeSettings(params.te_flag_overlap_frac, params.te_flag_evalue)
        nativeHelpers = ['te_search.py','prepare_protein_evidence.py'].collect { file("${params.pipeline_root}/py_scripts/${it}",checkIfExists:true) }
        BUILD_TE_REPEAT_DB(repeatRequests.map { aid,repeats,library -> tuple(aid,library) },nativeHelpers)
        SEARCH_TE_PROTEINS(annotations.map { aid,genome,gff,proteins,cds,validation,identity -> tuple(aid,proteins) }
            .join(BUILD_TE_REPEAT_DB.out.database),nativeHelpers)
        FLAG_TE_GENES(annotations.map { aid,genome,gff,proteins,cds,validation,identity -> tuple(aid,gff,identity) }
            .join(PREPARE_MODEL_METRICS.out.metrics).join(repeatRequests).join(SEARCH_TE_PROTEINS.out.search)
            .map { aid,gff,identity,metrics,regions,repeats,library,search -> tuple(aid,gff,metrics,regions,identity,repeats,search) },
            ['te_evidence.py','prepare_protein_evidence.py'].collect { file("${params.pipeline_root}/py_scripts/${it}",checkIfExists:true) })
        teReports = FLAG_TE_GENES.out.te_gene_report.join(FLAG_TE_GENES.out.status).toList()
    }
    optional = aedReports.toList().map { rows -> [rows] }
        .combine(junctionReports.toList().map { rows -> [rows] }).combine(teReports.map { rows -> [rows] })
    REPORT_MODEL_ASSESSMENT(PREPARE_MODEL_METRICS.out.metrics.combine(optional).map { aid,metrics,regions,aeds,junctions,tes ->
        def aed = aeds.find { it[0]==aid }
        def junction = junctions.find { it[0]==aid }
        def te = tes.find { it[0]==aid }
        if (params.run_aed && !aed) error "Missing AED assessment for ${aid}"
        if (params.run_junction_qc && !junction) error "Missing junction assessment for ${aid}"
        if (params.run_te_flagging && !te) error "Missing TE assessment for ${aid}"
        tuple(aid,metrics,aed ? [aed[1],aed[2]] : [],junction ? [junction[1],junction[2]] : [],te ? [te[1],te[2]] : [])
    },['model_assessment.py','prepare_protein_evidence.py'].collect { file("${params.pipeline_root}/py_scripts/${it}",checkIfExists:true) })
    emit:
    status = REPORT_MODEL_ASSESSMENT.out.status
    transcripts = REPORT_MODEL_ASSESSMENT.out.transcripts
}
