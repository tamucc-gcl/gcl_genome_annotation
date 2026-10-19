nextflow.enable.dsl = 2
include { FUNCTIONAL_ANNOTATION } from './workflows/functionalAnnotation.nf'

process BUILD_FUNCTIONAL_TEST_DATABASE {
    label 'diamond'
    storeDir "${params.diamond_database}"
    input:
    path reference
    output:
    path('synthetic.dmnd'), emit: database
    script:
    """
    diamond makedb --in '${reference}' --db synthetic --threads ${task.cpus}
    """
}

workflow {
    def root = file(params.functional_inputs_root).toAbsolutePath().toString()
    def manifest = new groovy.json.JsonSlurperClassic().parse(new File("${root}/preflight/resolved_inputs.json"))
    java.nio.file.Files.createDirectories(file(params.diamond_database))
    BUILD_FUNCTIONAL_TEST_DATABASE(file("${root}/reference.fa"))
    if (params.run_functional_annotation) {
        def jobs = Channel.fromList(manifest.assemblies).combine(BUILD_FUNCTIONAL_TEST_DATABASE.out.database).map { row,database ->
            def settings = [:]+row.settings
            settings.run_diamond = settings.run_diamond && params.run_diamond
            tuple(row.assembly_id,settings,file("${root}/structural_annotation/identified/${row.assembly_id}"))
        }
        FUNCTIONAL_ANNOTATION(jobs,Channel.empty())
    }
}
