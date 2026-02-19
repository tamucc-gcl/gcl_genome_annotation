#!/usr/bin/env nextflow

/*
========================================================================================
    GENOME ANNOTATION PIPELINE
========================================================================================
    Repeat modeling, curation, and masking pipeline

    Usage:
        nextflow run main.nf \
            --sample_sheet samples.csv \
            --dfam_dir /path/to/dfam \
            --taxa_name Actinopterygii \
            --outdir results \
            -profile singularity
========================================================================================
*/

nextflow.enable.dsl = 2

// ── Parameter defaults ───────────────────────────────────────────────────────
params.sample_sheet   = null
params.outdir         = "${launchDir}/results"
params.publish_dir_mode = 'link'

// Dfam
params.dfam_dir       = null
params.taxa_name      = null

// Repeat modeller
params.repModel_rounds = 3

// Repeat curation
params.min_repeat_length = 100
params.cluster_pct       = 0.95   // set to 0 or null to skip clustering

// RepeatMasker
params.job_speed = 's'            // null, 's', 'q', 'qq'

// Breaker Params
params.braker_rounds     = 5
params.braker_ab_initio  = false
params.protein_database  = null   // path to OrthoDB proteins

// Structural QC
params.busco_lineage = 'actinopterygii_odb10'
params.busco_downloads = '/work/birdlab/databases/busco_datasets' 
params.oma_db_name = 'LUCA'   
params.omark_downloads = '/work/birdlab/databases/omark_datasets' 
params.ncbi_speciesID = null

// ── Parameter validation ─────────────────────────────────────────────────────
if (!params.sample_sheet) { error "Please provide --sample_sheet" }
if (!params.dfam_dir)     { error "Please provide --dfam_dir" }
if (!params.taxa_name)    { error "Please provide --taxa_name" }

/*
========================================================================================
    IMPORT FUNCTIONS
========================================================================================
*/
include { parseSampleSheet } from './functions/parse_sample_sheet.nf'

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { REPEAT_MODEL  } from './modules/repeatModel.nf'
include { CURATE_REPEATS } from './modules/curateRepeats.nf'
include { CLASSIFY_REPEATS } from './modules/classifyRepeats.nf'
include { CURATE_DFAM   } from './modules/curateDfam.nf'
include { REPEAT_MASK   } from './modules/repeatMask.nf'
include { BRAKER } from './modules/braker.nf'
include { DOWNLOAD_BUSCO_DB } from './modules/download_busco_db.nf'
include { BUSCO } from './modules/busco.nf'
include { DOWNLOAD_OMARK_DB } from './modules/download_omark_db.nf'
include { OMARK } from './modules/omark.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/
workflow {

    // Parse sample sheet → channel of [assembly_id, assembly_path]
    ch_assemblies = parseSampleSheet(params.sample_sheet)

    /*
    ========================================================================================
        Database Downloads
    ========================================================================================
    */

    // ── Step 0: Extract Dfam library (runs once) ─────────────────────────────
    ch_dfam_dir = Channel.fromPath(params.dfam_dir, checkIfExists: true)
    CURATE_DFAM(ch_dfam_dir, params.taxa_name)

    // ── Step 0: BUSCO (runs once) ─────────────────────────────
    DOWNLOAD_BUSCO_DB(params.busco_lineage)
    ch_busco_db = DOWNLOAD_BUSCO_DB.out.db

    // ── Step 0: OMArk (runs once) ─────────────────────────────
    DOWNLOAD_OMARK_DB()

    /*
    ========================================================================================
        Repeat Annotation
    ========================================================================================
    */

    // ── Step 1: De novo repeat modeling (per genome) ─────────────────────────
    REPEAT_MODEL(ch_assemblies)

    // ── Step 2: Curate de novo repeats (per genome) ──────────────────────────
    CURATE_REPEATS(REPEAT_MODEL.out.families_fa)

    // ── Step 3: Classify de novo repeats (per genome) ──────────────────────────
    CLASSIFY_REPEATS(CURATE_REPEATS.out.clustered_repeats)

    // ── Step 4: Mask genome with combined libraries (per genome) ─────────────
    // Join assembly path back with curated repeats on assembly_id
    ch_mask_input = ch_assemblies
        .join(CLASSIFY_REPEATS.out.curated_repeats)
        // produces: [assembly_id, assembly_path, curated_repeats]

    REPEAT_MASK(ch_mask_input, 
                CURATE_DFAM.out.dfam_lib.collect())

    /*
    ========================================================================================
        Gene Structural Annotation
    ========================================================================================
    */

    // ── Step 5: Run Breaker on SoftMasked genome (per genome) ─────────────
    ch_prot_db = params.protein_database
        ? Channel.fromPath(params.protein_database, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_PROT_DB")
    //ch_augustus_config = Channel.fromPath(params.augustus_config, checkIfExists: true)

    BRAKER(REPEAT_MASK.out.masked_fasta, 
            ch_prot_db.collect()) 
            //ch_augustus_config.collect())

    // ── Step 6: Run BUSCO on protein predictions (per genome) ─────────────
    BUSCO(
        BRAKER.out.braker_proteins,
        DOWNLOAD_BUSCO_DB.out.db)

    // ── Step 7: Run OMArk on protein predictions (per genome) ─────────────
    OMARK(
        BRAKER.out.braker_proteins,
        DOWNLOAD_OMARK_DB.out.omamer_db,
        DOWNLOAD_OMARK_DB.out.taxa_db)

}

// ── On-complete summary ──────────────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ============================================================
    Pipeline complete!
    Results : ${params.outdir}
    Duration: ${workflow.duration}
    Status  : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    ============================================================
    """.stripIndent()
}