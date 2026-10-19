nextflow.enable.dsl = 2
include { GENOME_TRACKS } from './workflows/genomeTracks.nf'

// A locally generated covariance model exercises native cmscan without downloads.
process BUILD_TRACK_TEST_RFAM {
    label 'download_rfam_db'
    storeDir "${params.rfam_database}"
    input:
    path stockholm
    output:
    path('rfam_db'), emit: database
    script:
    """
    set -euo pipefail
    mkdir rfam_db
    cmbuild rfam_db/Rfam.cm '${stockholm}'
    cmcalibrate --cpu ${task.cpus} --seed 42 rfam_db/Rfam.cm
    cmpress rfam_db/Rfam.cm
    printf 'CL00001 Test_RNA\\n' > rfam_db/Rfam.clanin
    """
}

workflow {
    def root = file(params.genome_tracks_test_root).toAbsolutePath().toString()
    java.nio.file.Files.createDirectories(file(params.rfam_database))
    BUILD_TRACK_TEST_RFAM(file("${root}/model.sto"))
    def rows = new groovy.json.JsonSlurperClassic().parse(new File("${root}/assemblies.json"))
    def jobs = Channel.fromList(rows).combine(BUILD_TRACK_TEST_RFAM.out.database).map { row,db ->
        def settings = row.settings.collectEntries { k,v -> [k,v] }
        if (params.genome_tracks_test_disabled) ['trnascan','infernal','telomere','centromere'].each { settings['run_'+it]=false }
        tuple(row.assembly_id,settings,file("${root}/${row.assembly_id}.fa"),file("${root}/${row.assembly_id}.json"))
    }
    def masked = Channel.fromList(rows).map { row -> tuple(row.assembly_id,file("${root}/${row.assembly_id}.masked.fa")) }
    GENOME_TRACKS(jobs,masked)
}
