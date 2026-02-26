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
params.min_contig_length = 50000

// Dfam
params.dfam_version = '39'
params.dfam_database  = null
params.taxa_name      = null

// Repeat modeller
params.repModel_rounds = 3

// Repeat curation
params.min_repeat_length = 100
params.cluster_pct       = 0.95   // set to 0 or null to skip clustering

// RepeatMasker
params.job_speed = 's'            // null, 's', 'q', 'qq'
params.n_rm_chunks = 8          // number of chunks to split genome into for parallel RepeatMasker runs. Adjust based on genome size and available resources.

// Infernal
params.n_infernal_chunks = 8 // number of chunks to split genome into for parallel Infernal runs. Adjust based on genome size and available resources.
params.rfam_database = null // path to Rfam database directory (contains Rfam.cm and Rfam.clanin files)
params.rfam_version = 'CURRENT' // only one version

// Braker Params
params.braker_rounds     = 5
params.braker_ab_initio  = false
params.protein_database  = null   // path to OrthoDB proteins

// Structural QC
params.busco_lineage = 'actinopterygii_odb10'
params.busco_database = '/work/birdlab/databases/busco' 
params.oma_db_name = 'LUCA'   //https://omabrowser.org/oma/current/ see here for names
params.omark_database = '/work/birdlab/databases/omark' 
params.ncbi_speciesID = null

// Functional Annotations
params.n_chunks = 8 // number of chunks to split proteins into for parallel processing in InterProScan, Diamond, eggNOG. Adjust based on size of proteome and available resources.

params.interproscan_database = null //'/work/birdlab/databases/interproscan_datasets'
params.interproscan_version = '5.77-108.0' //'5.77-108.0'

params.eggnog_database = null //'/work/birdlab/databases/eggnog_datasets'
params.eggnog_version = '5.0.2' //'5.0.2'
params.eggnog_evalue = 0.001
params.eggnog_sensemode = 'sensitive' //{fast,mid-sensitive,sensitive,more-sensitive,very-sensitive,ultra-sensitive}
params.eggnog_taxid = 'none' // set to specific NCBI taxid (e.g. 7898 for Actinopterygii) to restrict to that lineage, or 'none' for no restriction
params.eggnog_target_orthologs = 'all'
params.eggnog_goEvidence = "non-electronic"

params.diamond_database = null // path to diamond database (e.g. nr or custom db)
params.diamond_version = 'swissprot' //only swissprot download enabled. Will get most recent version available at time of pipeline run. To use custom diamond database, set --diamond_database to path of database and set --diamond_version to custom
params.diamond_sensitivity = 'sensitive' // {fast, sensitive, more-sensitive, very-sensitive, ultra-sensitive}
params.diamond_maxHits = 1
params.diamond_evalue = 1e-5
params.diamond_qcov = 70

params.kofamscan_database = null // path to kofamscan database directory (contains profiles and ko_list)
params.kofamscan_version = 'latest' // only one version

// ── Parameter validation ─────────────────────────────────────────────────────
if (!params.sample_sheet) { error "Please provide --sample_sheet" }
if (!params.dfam_database) { error "Please provide --dfam_database" }
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

include { FILTER_ASSEMBLY } from './modules/filterAssembly.nf'

// ── Database ─────────────────────────────
include { CURATE_DFAM; DOWNLOAD_DFAM_DB } from './modules/curateDfam.nf'

// ── Repeat annotation ─────────────────────────────
include { REPEAT_MODEL  } from './modules/repeatModel.nf'
include { CURATE_REPEATS } from './modules/curateRepeats.nf'
include { CLASSIFY_REPEATS } from './modules/classifyRepeats.nf'
//include { REPEAT_MASK   } from './modules/repeatMask.nf'
include { REPEAT_MASK_PARALLEL as REPEAT_MASK } from './workflows/repeatMaskParallel.nf'

// ── Structural annotation ─────────────────────────────
include { BRAKER } from './modules/braker.nf'
include { AGAT } from './modules/agat.nf'
include { BUSCO; DOWNLOAD_BUSCO_DB } from './modules/busco.nf'
include { OMARK; DOWNLOAD_OMARK_DB } from './modules/omark.nf'

// ── ncRNA annotation ─────────────────────────────
include { tRNASCAN } from './modules/trna-scan.nf'
//include { INFERNAL; DOWNLOAD_RFAM_DB } from './modules/infernal.nf'
include { DOWNLOAD_RFAM_DB } from './modules/infernal.nf'
include { INFERNAL_PARALLEL as INFERNAL } from './workflows/infernalParallel.nf'

// ── Functional annotation ─────────────────────────────
include { CHUNK_PROTEINS } from './modules/chunk_proteins.nf'
include { INTERPROSCAN; MERGE_INTERPROSCAN; DOWNLOAD_INTERPROSCAN_DB } from './modules/interproscan.nf'
include { EGGNOG; MERGE_EGGNOG; DOWNLOAD_EGGNOG_DB } from './modules/eggnog.nf'
include { DIAMOND; MERGE_DIAMOND; DOWNLOAD_DIAMOND_DB } from './modules/diamond.nf'
include { KOFAMSCAN; MERGE_KOFAMSCAN; DOWNLOAD_KOFAMSCAN_DB } from './modules/kofamscan.nf'
include { DEEPSIG; MERGE_DEEPSIG } from './modules/deepsig.nf'
include { MERGE_FUNCTIONAL_ANNOTATIONS } from './modules/merge_functional_annotations.nf'

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

    // ── Step 0: DFam (runs once) ─────────────────────────────
    DOWNLOAD_DFAM_DB(params.dfam_version)

    // ── Step 0: RFam (runs once) ─────────────────────────────
    DOWNLOAD_RFAM_DB(params.rfam_version)

    // ── Step 0: BUSCO (runs once) ─────────────────────────────
    DOWNLOAD_BUSCO_DB(params.busco_lineage)
    ch_busco_db = DOWNLOAD_BUSCO_DB.out.db

    // ── Step 0: OMArk (runs once) ─────────────────────────────
    DOWNLOAD_OMARK_DB(params.oma_db_name)

    // ── Step 0: Download InterProScan database (runs once) ─────────────────────────────
    DOWNLOAD_INTERPROSCAN_DB(params.interproscan_version)

    // ── Step 0: Download EggNOG database (runs once) ─────────────────────────────
    DOWNLOAD_EGGNOG_DB(params.eggnog_version)

    // ── Step 0: Download Diamond database (runs once) ─────────────────────────────
    DOWNLOAD_DIAMOND_DB(params.diamond_version)

    // ── Step 0: Download KoFamScan database (runs once) ─────────────────────────────
    DOWNLOAD_KOFAMSCAN_DB(params.kofamscan_version)

    /*
    ========================================================================================
        Assembly Filtering and Stats
    ========================================================================================
    */

    FILTER_ASSEMBLY(ch_assemblies)

    ch_assemblies_filtered = FILTER_ASSEMBLY.out.filtered_assembly

    /*
    ========================================================================================
        Repeat Annotation
    ========================================================================================
    */

    // ── Step 0: Extract Dfam library (runs once) ─────────────────────────────
    CURATE_DFAM(DOWNLOAD_DFAM_DB.out.dfam_db_dir, 
                params.taxa_name)

    // ── Step 1: De novo repeat modeling (per genome) ─────────────────────────
    REPEAT_MODEL(ch_assemblies_filtered, 
                CURATE_DFAM.out.dfam_lib.collect())

    // ── Step 2: Curate de novo repeats (per genome) ──────────────────────────
    CURATE_REPEATS(REPEAT_MODEL.out.families_fa)

    // ── Step 3: Classify de novo repeats (per genome) ──────────────────────────
    CLASSIFY_REPEATS(CURATE_REPEATS.out.clustered_repeats)

    // ── Step 4: Mask genome with combined libraries (per genome) ─────────────
    // Join assembly path back with curated repeats on assembly_id
    ch_mask_input = ch_assemblies_filtered
        .join(CLASSIFY_REPEATS.out.curated_repeats)
        // produces: [assembly_id, assembly_path, curated_repeats]

    REPEAT_MASK(ch_mask_input, 
                CURATE_DFAM.out.dfam_lib.collect())

    /*
    ========================================================================================
        Structural Annotation
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

    // ── Step 6: Get descriptive stats for annotations (per genome) ─────────────
    ch_agat_inputs = ch_assemblies_filtered
        .join(BRAKER.out.braker_gff3, by: 0)
        .map { id, assembly_path, gff_path -> tuple(id, assembly_path, gff_path) }

    AGAT(ch_agat_inputs)

    // ── Step 7: Run BUSCO on protein predictions (per genome) ─────────────
    BUSCO(
        BRAKER.out.braker_proteins,
        DOWNLOAD_BUSCO_DB.out.db)

    // ── Step 8: Run OMArk on protein predictions (per genome) ─────────────
    OMARK(
        BRAKER.out.braker_proteins,
        DOWNLOAD_OMARK_DB.out.omamer_db,
        DOWNLOAD_OMARK_DB.out.taxa_db)

    /*
    ========================================================================================
        ncRNA Annotation
    ========================================================================================
    */

    // ── Step 9: Run tRNA-Scan-SE on SoftMasked genome (per genome) ─────────────
    tRNASCAN(REPEAT_MASK.out.masked_fasta)

    // ─ Step 10: Run Infernal on SoftMasked genome (per genome) ─────────────
    INFERNAL(REPEAT_MASK.out.masked_fasta, 
            DOWNLOAD_RFAM_DB.out.rfam_db_dir)

    /*
    ========================================================================================
        Functional Annotation
    ========================================================================================
    */

    // ── Step 11: Split proteins into chunks ─────────────
    CHUNK_PROTEINS(
        BRAKER.out.braker_proteins)

    // ch_chunks: (assembly_id, chunk_id, chunk_file)
    ch_chunks = CHUNK_PROTEINS.out.chunks
    .flatMap { assembly_id, f_or_list ->

        // Normalize: if it's already a list/collection, use it; otherwise wrap it
        def files = (f_or_list instanceof Collection) ? f_or_list : [f_or_list]

        files
            .sort { it.name }
            .collect { f ->
                def m = (f.name =~ /part_(\d+)/)
                def chunk_id = m ? "part_${m[0][1]}" : f.baseName
                tuple(assembly_id, chunk_id, f)
            }
    }
    // ── Step 12: Run InterProScan on protein chunks ─────────────
    INTERPROSCAN(
        ch_chunks,
        DOWNLOAD_INTERPROSCAN_DB.out.db)

    // group per assembly
    ch_ipr_tsv_grouped = INTERPROSCAN.out.interproscan_tsv
        .map { assembly_id, chunk_id, tsv -> tuple(assembly_id, tsv) }
        .groupTuple(size: params.n_chunks)

    MERGE_INTERPROSCAN(ch_ipr_tsv_grouped)

    // ── Step 13: Run eggNOG on protein chunks ─────────────
    EGGNOG(
        ch_chunks,
        DOWNLOAD_EGGNOG_DB.out.eggnog_db_dir)

    // group per assembly
    ch_eggnog_annotations_grouped = EGGNOG.out.eggnog_annotations
        .map { assembly_id, chunk_id, tsv -> tuple(assembly_id, tsv) }
        .groupTuple(size: params.n_chunks)

    MERGE_EGGNOG(ch_eggnog_annotations_grouped)

    // ── Step 14: Run Diamond on protein chunks ─────────────
    DIAMOND(
        ch_chunks,
        DOWNLOAD_DIAMOND_DB.out.diamond_db)

    // group per assembly
    ch_diamond_annotations_grouped = DIAMOND.out.diamond_annotations
        .map { assembly_id, chunk_id, tsv -> tuple(assembly_id, tsv) }
        .groupTuple(size: params.n_chunks)

    MERGE_DIAMOND(ch_diamond_annotations_grouped)

    // ── Step 15: Run KoFamScan on protein chunks ─────────────
    KOFAMSCAN(
        ch_chunks,
        DOWNLOAD_KOFAMSCAN_DB.out.kofamscan_db_dir)

    ch_kofamscan_grouped = KOFAMSCAN.out.kofamscan_annotations
        .map { assembly_id, chunk_id, tsv -> tuple(assembly_id, tsv) }
        .groupTuple(size: params.n_chunks)

    MERGE_KOFAMSCAN(ch_kofamscan_grouped)

    // ── Step 16: Run DeepSIG on protein chunks ─────────────
    DEEPSIG(ch_chunks)

    ch_deepsig_grouped = DEEPSIG.out.deepsig_annotations
        .map { assembly_id, chunk_id, gff -> tuple(assembly_id, gff) }
        .groupTuple(size: params.n_chunks)

    MERGE_DEEPSIG(ch_deepsig_grouped)

    // ── Step 17: Merge all functional annotations into single file per genome ─────────────
    ch_func_merged = MERGE_INTERPROSCAN.out.merged_tsv
        .join(MERGE_EGGNOG.out.merged)
        .join(MERGE_DIAMOND.out.merged)
        .join(MERGE_KOFAMSCAN.out.merged)
        .join(MERGE_DEEPSIG.out.merged_gff3)
        // produces: [assembly_id, ipr_tsv, eggnog_tsv, diamond_tsv, kofam_tsv, deepsig_gff3]

    MERGE_FUNCTIONAL_ANNOTATIONS(ch_func_merged)
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