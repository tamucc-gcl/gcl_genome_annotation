include { tRNASCAN } from '../modules/trna-scan.nf'
include { DOWNLOAD_RFAM_DB } from '../modules/infernal.nf'
include { INFERNAL_PARALLEL } from './infernalParallel.nf'
include { TELOMERE_ANNOTATION } from '../modules/telomere.nf'
include { CENTROMERE_ANNOTATION; CENTROMERE_TO_GFF3 } from '../modules/centromere.nf'
include { VALIDATE_GENOME_TRACK; HIGH_CONFIDENCE_TRNA; FILTER_GENOME_TRACK; REPORT_GENOME_TRACKS } from '../modules/genomeTracks.nf'

workflow GENOME_TRACKS {
    take:
    jobs // aid, resolved settings, filtered unmasked FASTA, identity dictionary
    masked // aid, softmasked FASTA
    main:
    if ((params.n_infernal_chunks as int) < 1) error 'n_infernal_chunks must be positive'
    if ((params.telomere_window as int) < 1 || (params.telomere_min_count as int) < 1 || (params.telomere_end_window as int) < 0)
        error 'Invalid telomere window/count/end distance'
    if (!Double.isFinite(params.telomere_filter_min_density as double) || (params.telomere_filter_min_density as double)<0 || (params.telomere_filter_min_density as double)>1 || (params.telomere_filter_max_end_distance as int)<0)
        error 'Telomere filter needs density 0..1 and nonnegative end distance'
    if (!Double.isFinite(params.centromere_filter_min_identity as double) || (params.centromere_filter_min_identity as double)<0 || (params.centromere_filter_min_identity as double)>1)
        error 'Centromere filter identity must be 0..1'
    [params.trna_filter_cmscore,params.trna_filter_ssscore,params.trna_filter_isoscore,params.trna_filter_tertiary_start,params.trna_filter_tertiary_max].each { value ->
        if (!Double.isFinite(value as double) || (value as double)< -1) error 'Invalid native tRNA confidence cutoff'
    }
    helpers = ['genome_tracks.py','reference_repeats.py'].collect { file("${projectDir}/py_scripts/${it}",checkIfExists:true) }
    reportHelpers = helpers + [file("${projectDir}/py_scripts/report_genome_tracks.py",checkIfExists:true)]
    tRNASCAN(jobs.filter { it[1].run_trnascan }.map { aid,settings,genome,identity -> tuple(aid,genome) })
    HIGH_CONFIDENCE_TRNA(tRNASCAN.out.confidence_inputs.filter { params.trna_filter })
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
    policies = [
        telomere:[enabled:params.telomere_filter,min_density:params.telomere_filter_min_density,max_end_distance:params.telomere_filter_max_end_distance,
                  raw_min_count:params.telomere_min_count,window:params.telomere_window],
        trnascan:[enabled:params.trna_filter,method:'EukHighConfidenceFilter',cmscore:params.trna_filter_cmscore,ssscore:params.trna_filter_ssscore,
                  isoscore:params.trna_filter_isoscore,tertiary_start:params.trna_filter_tertiary_start,tertiary_max:params.trna_filter_tertiary_max],
        infernal:[enabled:params.ncrna_filter,remove_overlaps:params.ncrna_filter_remove_overlaps,method:'native_inclusion_and_overlap'],
        centromere:[enabled:params.centromere_filter,min_identity:params.centromere_filter_min_identity,high_confidence_hor:params.centromere_filter_high_confidence_hor]
    ]
    auxiliary = TELOMERE_ANNOTATION.out.telomere_tsv.map { aid,path -> tuple(aid,'telomere',path) }
        .mix(HIGH_CONFIDENCE_TRNA.out.classifications.map { aid,folder -> tuple(aid,'trnascan',file("${folder}/${aid}.out",checkIfExists:true)) })
        .toList().map { rows -> [rows] }
    filtering = VALIDATE_GENOME_TRACK.out.validated.combine(jobs.map { aid,settings,genome,identity -> tuple(aid,identity,settings.telomere_motif) },by:0)
        .combine(auxiliary).map { aid,track,gff,status,identity,motif,rows ->
            def policy=[:]+policies[track]
            if (track=='telomere') policy.motif=motif
            tuple(aid,track,gff,status,identity,policy,policy.enabled ? rows.findAll { it[0]==aid && it[1]==track }.collect { it[2] } : [])
        }
    FILTER_GENOME_TRACK(filtering,[file("${projectDir}/py_scripts/filter_genome_tracks.py",checkIfExists:true)]+helpers)
    completed = FILTER_GENOME_TRACK.out.filtered.map { aid,track,gff,status,decisions -> tuple(aid,status) }.toList()
    REPORT_GENOME_TRACKS(jobs.combine(completed.map { rows -> [rows] }).map { aid,settings,genome,identity,rows ->
        def policy = settings.findAll { key,value -> key in ['run_trnascan','run_infernal','run_telomere','run_centromere','telomere_motif'] }
        policy += [rfam_version:params.rfam_version, telomere_window:params.telomere_window,
                   telomere_min_count:params.telomere_min_count,telomere_end_window:params.telomere_end_window]
        tuple(aid,policy,identity,rows.findAll { it[0] == aid }.collect { it[1] }.sort { it.name })
    },reportHelpers)
    emit:
    status = REPORT_GENOME_TRACKS.out.status
    gff3 = FILTER_GENOME_TRACK.out.filtered.map { aid,track,gff,status,decisions -> tuple(aid,track,gff) }
    raw_gff3 = VALIDATE_GENOME_TRACK.out.validated.map { aid,track,gff,status -> tuple(aid,track,gff) }
}
