include { PREPARE_FUNCTIONAL_INPUTS; EMPTY_FUNCTIONAL_RESULTS; REPORT_FUNCTIONAL_ANNOTATION } from '../modules/functionalAnnotation.nf'
include { CHUNK_PROTEINS } from '../modules/chunk_proteins.nf'
include { INTERPROSCAN; MERGE_INTERPROSCAN; DOWNLOAD_INTERPROSCAN_DB } from '../modules/interproscan.nf'
include { EGGNOG; MERGE_EGGNOG; DOWNLOAD_EGGNOG_DB } from '../modules/eggnog.nf'
include { DIAMOND; MERGE_DIAMOND; DOWNLOAD_DIAMOND_DB } from '../modules/diamond.nf'
include { KOFAMSCAN; MERGE_KOFAMSCAN; DOWNLOAD_KOFAMSCAN_DB } from '../modules/kofamscan.nf'
include { DEEPSIG; MERGE_DEEPSIG } from '../modules/deepsig.nf'
include { MERGE_FUNCTIONAL_ANNOTATIONS } from '../modules/merge_functional_annotations.nf'
include { ANNOTATE_GFF3 } from '../modules/annotate_gff3.nf'

workflow FUNCTIONAL_ANNOTATION {
    take:
    jobs // aid, resolved per-assembly settings, identified bundle
    omarkClasses // optional aid, native selected-representative .ump
    main:
    if ((params.n_chunks as int) < 1) error 'n_chunks must be positive'
    helpers = ['functional_inputs.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${projectDir}/py_scripts/${it}",checkIfExists:true)
    }
    PREPARE_FUNCTIONAL_INPUTS(jobs.map { aid,settings,folder -> tuple(aid,
        file("${folder}/${aid}_proteins.fa"),file("${folder}/${aid}.gff3"),
        file("${folder}/transcript_assessment.tsv"),file("${folder}/annotation_status.json"))
    },helpers)
    ready = jobs.join(PREPARE_FUNCTIONAL_INPUTS.out.prepared).map { aid,settings,folder,proteins,eligibility,validation ->
        def count = new groovy.json.JsonSlurperClassic().parse(new File(validation.toString())).eligible
        tuple(aid,settings,folder,proteins,eligibility,count)
    }
    CHUNK_PROTEINS(ready.filter { row -> row[5] > 0 && ['interproscan','eggnog','diamond','kofamscan','deepsig'].any { row[1]['run_'+it] } }
        .map { aid,settings,folder,proteins,eligibility,count -> tuple(aid,proteins) })
    chunks = CHUNK_PROTEINS.out.chunks.flatMap { aid,paths ->
        (paths instanceof List ? paths : [paths]).sort { it.name }.collect { path -> tuple(aid,path.baseName,path) }
    // Multiple chunks share an assembly key. join consumes a settings row once;
    // combine(by:0) attaches that row to every chunk for the assembly.
    }.combine(jobs.map { aid,settings,folder -> tuple(aid,settings) }, by: 0)

    // Shared sources are requested only by assemblies actually using each tool.
    // Existing storeDir outputs are reused even in a fresh project without -resume.
    iprRequests = chunks.filter { it[3].run_interproscan }.map { aid,cid,protein,settings -> tuple(aid,cid,protein) }
    iprVersions = iprRequests.map {
        if (!params.interproscan_database) error 'InterProScan requires --interproscan_database'
        java.nio.file.Files.createDirectories(file(params.interproscan_database))
        params.interproscan_version.toString()
    }.unique()
    DOWNLOAD_INTERPROSCAN_DB(iprVersions)
    INTERPROSCAN(iprRequests,DOWNLOAD_INTERPROSCAN_DB.out.db.first().map { path -> tuple(path.toAbsolutePath().toString(),path) })
    MERGE_INTERPROSCAN(INTERPROSCAN.out.interproscan_tsv.map { aid,cid,path -> tuple(aid,path) }.groupTuple())

    eggRequests = chunks.filter { it[3].run_eggnog }.map { aid,cid,protein,settings ->
        def scope = settings.eggnog_taxid ?: 'none'
        if (!(scope.toString() ==~ /(?:none|auto|[0-9]+(?:,[0-9]+)*)/)) error "Invalid eggNOG scope for ${aid}"
        tuple(aid,cid,protein,scope.toString())
    }
    eggVersions = eggRequests.map {
        if (!params.eggnog_database) error 'eggNOG requires --eggnog_database'
        java.nio.file.Files.createDirectories(file(params.eggnog_database))
        params.eggnog_version.toString()
    }.unique()
    DOWNLOAD_EGGNOG_DB(eggVersions)
    EGGNOG(eggRequests,DOWNLOAD_EGGNOG_DB.out.eggnog_db_dir.first())
    MERGE_EGGNOG(EGGNOG.out.eggnog_annotations.map { aid,cid,path -> tuple(aid,path) }.groupTuple())

    diamondRequests = chunks.filter { it[3].run_diamond }.map { aid,cid,protein,settings -> tuple(aid,cid,protein) }
    diamondVersions = diamondRequests.map {
        if (!params.diamond_database) error 'DIAMOND requires --diamond_database (shared store directory)'
        java.nio.file.Files.createDirectories(file(params.diamond_database))
        if (params.diamond_version != 'swissprot' && !file("${params.diamond_database}/${params.diamond_version}.dmnd").exists())
            error 'Custom DIAMOND index is missing from diamond_database'
        params.diamond_version.toString()
    }.unique()
    DOWNLOAD_DIAMOND_DB(diamondVersions)
    DIAMOND(diamondRequests,DOWNLOAD_DIAMOND_DB.out.diamond_db.first())
    MERGE_DIAMOND(DIAMOND.out.diamond_annotations.map { aid,cid,path -> tuple(aid,path) }.groupTuple())

    kofamRequests = chunks.filter { it[3].run_kofamscan }.map { aid,cid,protein,settings -> tuple(aid,cid,protein) }
    kofamVersions = kofamRequests.map {
        if (!params.kofamscan_database) error 'KofamScan requires --kofamscan_database'
        java.nio.file.Files.createDirectories(file(params.kofamscan_database))
        params.kofamscan_version.toString()
    }.unique()
    DOWNLOAD_KOFAMSCAN_DB(kofamVersions)
    KOFAMSCAN(kofamRequests,DOWNLOAD_KOFAMSCAN_DB.out.kofamscan_db_dir.first())
    MERGE_KOFAMSCAN(KOFAMSCAN.out.kofamscan_annotations.map { aid,cid,path -> tuple(aid,path) }.groupTuple())

    DEEPSIG(chunks.filter { it[3].run_deepsig }.map { aid,cid,protein,settings -> tuple(aid,cid,protein) })
    MERGE_DEEPSIG(DEEPSIG.out.deepsig_annotations.map { aid,cid,path -> tuple(aid,path) }.groupTuple())
    EMPTY_FUNCTIONAL_RESULTS(jobs.map { aid,settings,folder -> aid })
    results = MERGE_INTERPROSCAN.out.merged_tsv.map { aid,path -> tuple(aid,'interproscan',path) }
        .mix(MERGE_EGGNOG.out.merged.map { aid,path -> tuple(aid,'eggnog',path) })
        .mix(MERGE_DIAMOND.out.merged.map { aid,path -> tuple(aid,'diamond',path) })
        .mix(MERGE_KOFAMSCAN.out.merged.map { aid,path -> tuple(aid,'kofamscan',path) })
        .mix(MERGE_DEEPSIG.out.merged_gff3.map { aid,path -> tuple(aid,'deepsig',path) }).toList().map { [it] }
    optionalQc = omarkClasses.toList().map { [it] }
    mergeJobs = ready.join(EMPTY_FUNCTIONAL_RESULTS.out.results).combine(results).combine(optionalQc).map {
        aid,settings,folder,proteins,eligibility,count,empty,rows,qc ->
        def paths = ['interproscan','eggnog','diamond','kofamscan','deepsig'].collect { tool ->
            def match = rows.find { it[0]==aid && it[1]==tool }
            if (settings['run_'+tool] && count>0 && !match) error "Missing ${tool} results for ${aid}"
            match ? match[2] : file("${empty}/${tool}.${tool=='deepsig' ? 'gff3' : 'tsv'}")
        }
        def omark = qc.find { it[0]==aid }
        def provenance = [:]+settings
        provenance.functional_parameters = [n_chunks:params.n_chunks,
            interproscan_version:params.interproscan_version,eggnog_version:params.eggnog_version,
            eggnog_evalue:params.eggnog_evalue,eggnog_sensemode:params.eggnog_sensemode,
            eggnog_target_orthologs:params.eggnog_target_orthologs,eggnog_goEvidence:params.eggnog_goEvidence,
            diamond_version:params.diamond_version,diamond_evalue:params.diamond_evalue,
            diamond_sensitivity:params.diamond_sensitivity,diamond_maxHits:params.diamond_maxHits,
            diamond_qcov:params.diamond_qcov,kofamscan_version:params.kofamscan_version]
        provenance.database_stores = [interproscan:params.interproscan_database,eggnog:params.eggnog_database,
            diamond:params.diamond_database,kofamscan:params.kofamscan_database]
        def encoded = groovy.json.JsonOutput.toJson(provenance).getBytes('UTF-8').encodeBase64().toString()
        tuple(aid,encoded,file("${folder}/transcript_assessment.tsv"),eligibility,*paths,omark ? omark[1] : [])
    }
    MERGE_FUNCTIONAL_ANNOTATIONS(mergeJobs,['merge_functional_annotations.py','functional_inputs.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${projectDir}/py_scripts/${it}",checkIfExists:true)
    })
    ANNOTATE_GFF3(jobs.join(MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv).map { aid,settings,folder,table ->
        tuple(aid,file("${folder}/${aid}.gff3"),[],[],[],[],[],table)
    },file("${projectDir}/py_scripts/annotate_gff3.py",checkIfExists:true))
    REPORT_FUNCTIONAL_ANNOTATION(jobs.join(ANNOTATE_GFF3.out.annotated_gff3).join(MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv)
        .join(MERGE_FUNCTIONAL_ANNOTATIONS.out.status).join(PREPARE_FUNCTIONAL_INPUTS.out.prepared).map {
            aid,settings,folder,gff,table,status,proteins,eligibility,validation ->
            tuple(aid,file("${folder}/${aid}.gff3"),gff,table,status,eligibility)
        },['functional_report.py','functional_inputs.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
            file("${projectDir}/py_scripts/${it}",checkIfExists:true)
        })
    emit:
    annotation = ANNOTATE_GFF3.out.annotated_gff3
    table = MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv
    status = REPORT_FUNCTIONAL_ANNOTATION.out.status
    report = REPORT_FUNCTIONAL_ANNOTATION.out.report
}
