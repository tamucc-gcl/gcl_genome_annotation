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

// Containers
params.repeatmodeler_container = null  // e.g. '/path/to/repeatmodeler.sif'

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
include { CURATE_DFAM   } from './modules/curateDfam.nf'
include { REPEAT_MASK   } from './modules/repeatMask.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/
workflow {

    // Parse sample sheet → channel of [assembly_id, assembly_path]
    ch_assemblies = parseSampleSheet(params.sample_sheet)

    // ── Step 0: Extract Dfam library (runs once) ─────────────────────────────
    ch_dfam_dir = Channel.fromPath(params.dfam_dir, checkIfExists: true)
    CURATE_DFAM(ch_dfam_dir, params.taxa_name)

    // ── Step 1: De novo repeat modeling (per genome) ─────────────────────────
    REPEAT_MODEL(ch_assemblies)

    // ── Step 2: Curate de novo repeats (per genome) ──────────────────────────
    CURATE_REPEATS(REPEAT_MODEL.out.families_fa)

    // ── Step 3: Mask genome with combined libraries (per genome) ─────────────
    // Join assembly path back with curated repeats on assembly_id
    ch_mask_input = ch_assemblies
        .join(CURATE_REPEATS.out.curated_repeats)
        // produces: [assembly_id, assembly_path, curated_repeats]

    REPEAT_MASK(ch_mask_input, CURATE_DFAM.out.dfam_lib)
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