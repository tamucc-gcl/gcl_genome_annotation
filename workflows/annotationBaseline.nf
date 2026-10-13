// Reuse the legacy scientific modules with assembly-specific settings.
include { PREPARE_ANNOTATION_BASELINE; REPORT_ANNOTATION_BASELINE } from '../modules/annotationBaseline.nf'
include { AGAT } from '../modules/agat.nf'
include { BUSCO; DOWNLOAD_BUSCO_DB } from '../modules/busco.nf'
include { OMARK; DOWNLOAD_OMARK_DB } from '../modules/omark.nf'

workflow ANNOTATION_BASELINE {
    take:
    models // [aid, settings, genome, gff, proteins, validation]
    helpers

    main:
    PREPARE_ANNOTATION_BASELINE(models.map { aid, settings, genome, gff, proteins, validation ->
        tuple(aid, settings, genome, gff, proteins, validation)
    }, helpers)
    def ready = models.join(PREPARE_ANNOTATION_BASELINE.out.selection)
    AGAT(ready.map { aid, settings, genome, gff, proteins, validation, selection -> tuple(aid, genome, gff) })
    def buscoModels = ready.filter { it[1].run_busco }
    def lineages = buscoModels.map { it[1].busco_lineage.toString() }.unique()
    def lineageRequests = lineages.map { lineage ->
        java.nio.file.Files.createDirectories(file(params.busco_database))
        lineage
    }
    DOWNLOAD_BUSCO_DB(lineageRequests)
    def buscoRequests = buscoModels.map { aid, settings, genome, gff, proteins, validation, selection ->
        tuple(settings.busco_lineage.toString(), aid, file("${selection}/${aid}.fa", checkIfExists: true))
    }.combine(DOWNLOAD_BUSCO_DB.out.db, by: 0).map { lineage, aid, proteins, database ->
        tuple(aid, lineage, proteins, database)
    }
    BUSCO(buscoRequests)
    def omarkModels = ready.filter { it[1].run_omark }
    DOWNLOAD_OMARK_DB(omarkModels.map { params.oma_db_name.toString() }.unique().map { name ->
        java.nio.file.Files.createDirectories(file(params.omark_database))
        name
    })
    def omarkRequests = omarkModels.map { aid, settings, genome, gff, proteins, validation, selection ->
        tuple(aid, settings.taxid.toString(), file("${selection}/${aid}.fa", checkIfExists: true))
    }.combine(DOWNLOAD_OMARK_DB.out.databases).map { aid, taxid, proteins, name, db, taxonomy ->
        tuple(aid, taxid, proteins, db, taxonomy)
    }
    OMARK(omarkRequests)
    // Sort each collection once, carrying its assembly IDs with its paths.
    // Empty optional channels become bound empty lists, never an inner join.
    def agatRows = AGAT.out.annotation_stats.toList().map { rows -> rows.sort { it[0] } }
    def buscoRows = BUSCO.out.busco.toList().map { rows -> rows.sort { it[0] } }
    def omarkRows = OMARK.out.omark_detailed.toList().map { rows -> rows.sort { it[0] } }
    REPORT_ANNOTATION_BASELINE(agatRows.map { it.collect { row -> row[0] } },
        PREPARE_ANNOTATION_BASELINE.out.selection.map { aid, folder -> file("${folder}/selection.json") }.toList(),
        agatRows.map { it.collect { row -> row[1] } },
        buscoRows.map { it.collect { row -> row[0] } }, buscoRows.map { it.collect { row -> row[1] } },
        omarkRows.map { it.collect { row -> row[0] } }, omarkRows.map { it.collect { row -> row[1] } }, helpers)

    emit:
    status = REPORT_ANNOTATION_BASELINE.out.status
    report = REPORT_ANNOTATION_BASELINE.out.report
}
