include { FILTER_GENE_MODELS } from '../modules/filterGeneModels.nf'

workflow MODEL_SELECTION {
    take:
    annotations // aid, GFF, proteins, CDS, final validation, reference identity
    assessments // aid, transcript assessment, assessment status
    main:
    selectionSettings = groovy.json.JsonOutput.toJson([
        structural_action:params.filter_structural_action,
        min_cds_length:params.filter_min_cds_length,
        min_gene_length:params.filter_min_gene_length,
        min_intron_length:params.filter_min_intron_length,
        aed_action:params.filter_aed_action, aed_source:params.filter_aed_source, max_aed:params.filter_max_aed,
        junction_action:params.filter_junction_action, min_junction_fraction:params.filter_min_junction_fraction,
        te_action:params.filter_te_action
    ]).bytes.encodeBase64().toString()
    FILTER_GENE_MODELS(annotations.join(assessments),selectionSettings,
        ['model_selection.py','model_metrics.py','annotation_models.py','prepare_protein_evidence.py','reference_repeats.py'].collect {
            file("${projectDir}/py_scripts/${it}",checkIfExists:true)
        })
    emit:
    annotation = FILTER_GENE_MODELS.out.bundle
}
