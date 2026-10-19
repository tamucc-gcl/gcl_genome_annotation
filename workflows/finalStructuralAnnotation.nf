include { RENAME_GENE_MODELS } from '../modules/renameGeneModels.nf'
include { ANNOTATION_BASELINE as SELECTED_ANNOTATION_QC } from './annotationBaseline.nf'
include { STRUCTURAL_QC; STRUCTURAL_QC_INDEX } from '../modules/structural_qc.nf'

workflow FINAL_STRUCTURAL_ANNOTATION {
    take:
    metadata // aid, QC settings, genome, identity
    selected // aid, selected bundle
    assessments // aid, original unfiltered assessment, status
    beforeQc // aid, path; original QC files and representative selection metadata
    main:
    if (params.final_id_prefix && !(params.final_id_prefix ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/))
        error 'final_id_prefix must use letters, digits, underscore, dot or hyphen'
    RENAME_GENE_MODELS(selected.join(metadata).map { aid,folder,settings,genome,identity -> tuple(aid,folder,identity) },
        params.final_id_prefix ?: '',
        ['final_identifiers.py','model_selection.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
            file("${projectDir}/py_scripts/${it}",checkIfExists:true)
        })
    SELECTED_ANNOTATION_QC(metadata.join(RENAME_GENE_MODELS.out.bundle).map { aid,settings,genome,identity,folder ->
        tuple(aid,settings,genome,file("${folder}/${aid}.gff3"),file("${folder}/${aid}_proteins.fa"),file("${folder}/annotation_status.json"))
    },['identified_qc.py','annotation_baseline.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${projectDir}/py_scripts/${it}",checkIfExists:true)
    })
    afterQc = SELECTED_ANNOTATION_QC.out.agat.mix(SELECTED_ANNOTATION_QC.out.busco)
        .mix(SELECTED_ANNOTATION_QC.out.omark).mix(SELECTED_ANNOTATION_QC.out.omark_classes)
        .mix(SELECTED_ANNOTATION_QC.out.selections.map { aid,folder -> tuple(aid,file("${folder}/selection.json")) })
        .mix(SELECTED_ANNOTATION_QC.out.selections.map { aid,folder -> tuple(aid,file("${folder}/representatives.tsv")) })
    qcRows = beforeQc.toList().map { rows -> [rows] }.combine(afterQc.toList().map { rows -> [rows] })
    STRUCTURAL_QC(RENAME_GENE_MODELS.out.bundle.join(selected).join(assessments).combine(qcRows).map {
        aid,identified,selection,assessment,status,before,after ->
        tuple(aid,identified,selection,assessment,status,before.findAll { it[0]==aid }.collect { it[1] },after.findAll { it[0]==aid }.collect { it[1] })
    },['structural_report.py','collect_annotation_metrics.py','model_selection.py','model_metrics.py','evidence_aed.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
        file("${projectDir}/py_scripts/${it}",checkIfExists:true)
    })
    reportRows = STRUCTURAL_QC.out.report.toList().map { rows -> rows.sort { it[0] } }
    STRUCTURAL_QC_INDEX(metadata.toList().map { it.collect { row -> row[0] }.sort() },reportRows.map { it.collect { row -> row[1] } },
        file("${projectDir}/py_scripts/structural_report_index.py",checkIfExists:true))
    emit:
    annotation = RENAME_GENE_MODELS.out.bundle
    omark_classes = SELECTED_ANNOTATION_QC.out.omark_classes
    report = STRUCTURAL_QC.out.report
}
