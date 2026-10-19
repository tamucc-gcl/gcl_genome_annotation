include { tRNASCAN } from '../modules/trna-scan.nf'
include { DOWNLOAD_RFAM_DB } from '../modules/infernal.nf'
include { INFERNAL_PARALLEL } from './infernalParallel.nf'
include { TELOMERE_ANNOTATION } from '../modules/telomere.nf'
include { CENTROMERE_ANNOTATION; CENTROMERE_TO_GFF3 } from '../modules/centromere.nf'
include { VALIDATE_GENOME_TRACK; REPORT_GENOME_TRACKS } from '../modules/genomeTracks.nf'

workflow GENOME_TRACKS {
    take:
    jobs // aid, resolved settings, filtered unmasked FASTA, identity dictionary
    masked // aid, softmasked FASTA
    main:
    if ((params.n_infernal_chunks as int) < 1) error 'n_infernal_chunks must be positive'
    if ((params.telomere_window as int) < 1 || (params.telomere_min_count as int) < 1 || (params.telomere_end_window as int) < 0)
        error 'Invalid telomere window/count/end distance'
    helpers = ['genome_tracks.py','reference_repeats.py'].collect { file("${projectDir}/py_scripts/${it}",checkIfExists:true) }
    reportHelpers = helpers + [file("${projectDir}/py_scripts/report_genome_tracks.py",checkIfExists:true)]
    tRNASCAN(jobs.filter { it[1].run_trnascan }.map { aid,settings,genome,identity -> tuple(aid,genome) })
    infernalRequests = jobs.filter { it[1].run_infernal }.join(masked).map { aid,settings,genome,identity,softmasked -> tuple(aid,softmasked) }
    rfamVersions = infernalRequests.map {
        if (!params.rfam_database) error 'Infernal requires --rfam_database (shared Rfam store directory)'
        java.nio.file.Files.createDirectories(file(params.rfam_database))
        params.rfam_version.toString()
    }.unique()
    DOWNLOAD_RFAM_DB(rfamVersions)
    INFERNAL_PARALLEL(infernalRequests,DOWNLOAD_RFAM_DB.out.rfam_db_dir.first())
    TELOMERE_ANNOTATION(jobs.filter { it[1].run_telomere }.map { aid,settings,genome,identity ->
        if (!(settings.telomere_motif ==~ /[ACGT]+/)) error "Invalid resolved telomere motif for ${aid}"
        tuple(aid,genome,settings.telomere_motif)
    },helpers)
    CENTROMERE_ANNOTATION(jobs.filter { it[1].run_centromere }.map { aid,settings,genome,identity -> tuple(aid,genome) })
    CENTROMERE_TO_GFF3(CENTROMERE_ANNOTATION.out.monomer_csv.join(CENTROMERE_ANNOTATION.out.hor_csv),helpers)
    raw = tRNASCAN.out.trnascan_gff3.map { aid,path -> tuple(aid,'trnascan',path) }
        .mix(INFERNAL_PARALLEL.out.infernal_gff3.map { aid,path -> tuple(aid,'infernal',path) })
        .mix(TELOMERE_ANNOTATION.out.telomere_gff3.map { aid,path -> tuple(aid,'telomere',path) })
        .mix(CENTROMERE_TO_GFF3.out.centromere_gff3.map { aid,path -> tuple(aid,'centromere',path) })
    // Several tracks share each assembly key: combine broadcasts its dictionary to all tracks.
    VALIDATE_GENOME_TRACK(raw.combine(jobs.map { aid,settings,genome,identity -> tuple(aid,identity) },by:0),helpers)
    completed = VALIDATE_GENOME_TRACK.out.validated.map { aid,track,gff,status -> tuple(aid,status) }.toList()
    REPORT_GENOME_TRACKS(jobs.combine(completed.map { rows -> [rows] }).map { aid,settings,genome,identity,rows ->
        def policy = settings.findAll { key,value -> key in ['run_trnascan','run_infernal','run_telomere','run_centromere','telomere_motif'] }
        policy += [rfam_version:params.rfam_version, telomere_window:params.telomere_window,
                   telomere_min_count:params.telomere_min_count,telomere_end_window:params.telomere_end_window]
        tuple(aid,policy,identity,rows.findAll { it[0] == aid }.collect { it[1] }.sort { it.name })
    },reportHelpers)
    emit:
    status = REPORT_GENOME_TRACKS.out.status
    gff3 = VALIDATE_GENOME_TRACK.out.validated.map { aid,track,gff,status -> tuple(aid,track,gff) }
}
