nextflow.enable.dsl = 2
include { ANNOTATION_REPORT } from './modules/annotation_report.nf'

workflow {
    def request=new groovy.json.JsonSlurperClassic().parse(new File(params.report_test_request))
    ANNOTATION_REPORT(request.releases.collect { file(it) },request.structural.collect { it[0] },request.structural.collect { file(it[1]) },
        request.stages.collect { [it[0],it[1]] },request.stages.collect { file(it[2]) },request.configuration,
        ['generate_annotation_report.py','unified_report.py','reference_repeats.py'].collect { file("${projectDir}/py_scripts/${it}",checkIfExists:true) })
}
