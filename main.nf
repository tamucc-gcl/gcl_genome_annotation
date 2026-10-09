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

// Parameter defaults and process resources are defined in nextflow.config.

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

include { COMBINE_REPEAT_LIBS } from './modules/combineRepeatLibs.nf'
include { FILTER_ASSEMBLY } from './modules/filterAssembly.nf'

// ── Database ─────────────────────────────
include { CURATE_DFAM; DOWNLOAD_DFAM_DB } from './modules/curateDfam.nf'

// ── Repeat annotation ─────────────────────────────
include { REPEAT_MODEL  } from './modules/repeatModel.nf'
include { CURATE_REPEATS } from './modules/curateRepeats.nf'
include { CLASSIFY_REPEATS } from './modules/classifyRepeats.nf'
include { RECLASSIFY_REPEATS } from './modules/reclassifyRepeats.nf'
include { EXTRACT_REPEAT_PEPS; INDEX_REPEAT_PEPS; FILTER_REPEAT_LIBRARY } from './modules/filterRepeatLibrary.nf'
//include { REPEAT_MASK   } from './modules/repeatMask.nf'
include { REPEAT_MASK_PARALLEL as REPEAT_MASK } from './workflows/repeatMaskParallel.nf'

// ── Telomere annotation ─────────────────────────────
include { TELOMERE_ANNOTATION } from './modules/telomere.nf'

// ── Centromere annotation ─────────────────────────────
include { CENTROMERE_ANNOTATION; CENTROMERE_TO_GFF3 } from './modules/centromere.nf'

// ── Structural annotation ─────────────────────────────
include { DOWNLOAD_ORTHODB } from './modules/downloadOrthoDB.nf'
include { PREPARE_PROTEIN_DB } from './modules/prepareProteinDb.nf'
include { BRAKER } from './modules/braker.nf'
include { COMPUTE_AED; AED_STUB } from './modules/computeAed.nf'
include { FILTER_GENE_MODELS } from './modules/filterGeneModels.nf'
include { FLAG_TE_GENES } from './modules/flagTeGenes.nf'
include { RENAME_GENE_MODELS } from './modules/renameGeneModels.nf'
include { AGAT } from './modules/agat.nf'
include { BUSCO; DOWNLOAD_BUSCO_DB } from './modules/busco.nf'
include { OMARK; DOWNLOAD_OMARK_DB } from './modules/omark.nf'
include { STRUCTURAL_QC } from './modules/structural_qc.nf'

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

// ── Annotation merging ─────────────────────────────
include { MERGE_FUNCTIONAL_ANNOTATIONS } from './modules/merge_functional_annotations.nf'
include { ANNOTATE_GFF3 } from './modules/annotate_gff3.nf'
include { EXTRACT_SEQUENCES } from './modules/extractSequences.nf'
include { ANNOTATION_REPORT } from './modules/annotation_report.nf'


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

    // ── Step 0: Prepare protein database for BRAKER ─────────────────────────
    //   Supports four configurations:
    //     --orthodb_clade Vertebrata --ref_proteome_accession GCF_xxx → combined
    //     --orthodb_clade Vertebrata                                  → OrthoDB only
    //     --ref_proteome_accession GCF_xxx                            → reference only
    //     (neither specified)                                         → no protein evidence
    if (params.orthodb_clade) {
        DOWNLOAD_ORTHODB(params.orthodb_clade)
        ch_orthodb = DOWNLOAD_ORTHODB.out.orthodb_proteins
    } else {
        ch_orthodb = Channel.fromPath("${projectDir}/assets/NO_PROT_DB")
    }

    PREPARE_PROTEIN_DB(DOWNLOAD_ORTHODB.out.orthodb_proteins)

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

    // ── Step 0: Extract repeat peptides (runs once) ─────────────────────────────
    EXTRACT_REPEAT_PEPS()
    INDEX_REPEAT_PEPS(EXTRACT_REPEAT_PEPS.out.repeat_peps.map { path ->
        tuple('repeat_peps_' + java.security.MessageDigest.getInstance('SHA-256').digest(path.bytes).encodeHex().toString(), path)
    })

    /*
    ========================================================================================
        Assembly Filtering and Stats
    ========================================================================================
    */

    def reference_helper = file("${projectDir}/py_scripts/reference_repeats.py", checkIfExists: true)
    FILTER_ASSEMBLY(ch_assemblies.map { aid, fa -> tuple(aid, fa, params.min_contig_length) }, reference_helper)
    ch_assemblies_filtered = FILTER_ASSEMBLY.out.bundle.map { aid, folder ->
        def state = new groovy.json.JsonSlurper().parse(new File("${folder}/${aid}.reference.json"))
        if (state.status != 'ready') error "Reference preparation failed for ${aid}: ${state.reason}"
        tuple(aid, file("${folder}/${aid}.filtered.fasta"))
    }

    /*
    ========================================================================================
        Repeat Annotation
    ========================================================================================
    */

    // ── Step 0: Extract Dfam library (runs once) ─────────────────────────────
    def legacy_taxon_key = 'taxon_' + java.security.MessageDigest.getInstance('SHA-256')
        .digest(params.taxa_name.toString().getBytes('UTF-8')).encodeHex().toString().substring(0, 20)
    CURATE_DFAM(DOWNLOAD_DFAM_DB.out.dfam_db_dir.map { folder -> tuple(legacy_taxon_key, params.taxa_name, folder) })
    ch_dfam_library = CURATE_DFAM.out.bundle.map { key, folder -> file("${folder}/dfam_subset.lib") }
    ch_dfam_source = CURATE_DFAM.out.bundle.map { key, folder -> file("${folder}/dfam_subset.source.json") }

    // ── Step 1: De novo repeat modeling (per genome) ─────────────────────────
    REPEAT_MODEL(ch_assemblies_filtered.combine(ch_dfam_source), params.repeat_seed)

    // ── Step 2: Curate de novo repeats (per genome) ──────────────────────────
    CURATE_REPEATS(REPEAT_MODEL.out.families_fa, reference_helper)

    // ── Step 3: Classify de novo repeats (per genome) ──────────────────────────
    CLASSIFY_REPEATS(CURATE_REPEATS.out.clustered_repeats)

    // ── Step 3b: Reclassify Unknown families with Terrier (per genome) ───────
    //   Uses deep learning to assign superfamily labels to families that
    //   RepeatClassifier couldn't classify via homology. Already-classified
    //   families pass through unchanged.
    RECLASSIFY_REPEATS(CLASSIFY_REPEATS.out.curated_repeats)

    // ── Step 4: Filter gene-derived sequences from repeat library ──────────────
    //   Uses two-database strategy: RepeatPeps (TE proteins) + SwissProt
    //   Repeats hitting RepeatPeps are confirmed TEs → kept
    //   Repeats hitting only SwissProt above threshold → gene-derived → removed
    FILTER_REPEAT_LIBRARY(
        RECLASSIFY_REPEATS.out.reclassified_lib,
        DOWNLOAD_DIAMOND_DB.out.diamond_db,
        EXTRACT_REPEAT_PEPS.out.repeat_peps,
        INDEX_REPEAT_PEPS.out.index.map { key, db, source -> db })

    // ── Step 5: Mask genome with combined libraries (per genome) ─────────────
    // Join assembly path back with filtered repeats on assembly_id
    COMBINE_REPEAT_LIBS(FILTER_REPEAT_LIBRARY.out.filtered_lib.combine(ch_dfam_library))
    REPEAT_MASK(ch_assemblies_filtered.join(COMBINE_REPEAT_LIBS.out.combined_lib), params.n_rm_chunks)

    /*
    ========================================================================================
        Telomere Annotation
    ========================================================================================
    */
    TELOMERE_ANNOTATION(ch_assemblies_filtered)

    /*
    ========================================================================================
        Centromere Annotation
    ========================================================================================
    */

    // ── Step: De novo centromere detection with centroAnno ────────────────────
    //   Runs on UNMASKED assembly — centromeric satellites are tandem repeats
    //   that RepeatMasker soft-masks, so centroAnno must see the original
    //   sequence to detect them via k-mer frequency analysis.
    //   No dependency on repeat pipeline — starts immediately.
    CENTROMERE_ANNOTATION(ch_assemblies_filtered)

    // ── Step: Convert centroAnno output to GFF3 ──────────────────────────────
    ch_centroanno_gff_input = CENTROMERE_ANNOTATION.out.monomer_csv
        .join(CENTROMERE_ANNOTATION.out.hor_csv, by: 0)
        // produces: [assembly_id, monomer_csv, hor_csv]

    CENTROMERE_TO_GFF3(ch_centroanno_gff_input)


    /*
    ========================================================================================
        Structural Annotation
    ========================================================================================
    */

    // ── Step 5a: Run Breaker on SoftMasked genome (per genome) ─────────────
    BRAKER(REPEAT_MASK.out.masked_fasta,
            PREPARE_PROTEIN_DB.out.combined_proteins.collect())

    // ── Compute AED (runs immediately after BRAKER) ──────────────
    // protein_db is the combined FASTA from PREPARE_PROTEIN_DB.
    // Only compute AED if a protein database was actually built.

    ch_braker_for_aed = BRAKER.out.braker_gff3
        .join(REPEAT_MASK.out.masked_fasta)
        .map { aid, gff, fa -> tuple(aid, fa, gff) }

    if (params.orthodb_clade || params.ref_proteome_accession) {
        COMPUTE_AED(ch_braker_for_aed,
                PREPARE_PROTEIN_DB.out.combined_proteins.collect())
        ch_aed_scores = COMPUTE_AED.out.aed_scores
    } else {
        // No protein DB → skip AED, emit header-only stub via /dev/null trick
        AED_STUB(BRAKER.out.braker_gff3)
        ch_aed_scores = AED_STUB.out.aed_scores
    }

    // ── Step 5b: Run gene filtering (per genome) ─────────────
    FILTER_GENE_MODELS(
        BRAKER.out.braker_gff3
            .join(BRAKER.out.braker_proteins)
            .join(ch_aed_scores)
            // → (assembly_id, gff3, proteins, aed_scores)
        )


    // ── Step 5c: Flag TE-derived gene models ───────────────────────────────
    ch_flag_input = FILTER_GENE_MODELS.out.filtered_gff3
        .join(FILTER_GENE_MODELS.out.filtered_proteins)
        .join(REPEAT_MASK.out.repeat_gff)
        .join(RECLASSIFY_REPEATS.out.reclassified_lib)
        // produces: [assembly_id, braker_gff3, braker_proteins, repeat_gff]

    FLAG_TE_GENES(ch_flag_input)
    //   ↑ Use the PRE-filter repeat library here — we want all TE sequences
    //     for homology detection, including ones removed from masking library

    // ── Rename retained protein-coding gene models ──────────────────────────────
    // This is the canonical ID boundary.
    // After this point, all downstream protein-coding IDs are final IDs.

    ch_rename_input = FLAG_TE_GENES.out.flagged_gff3
        .join(FILTER_GENE_MODELS.out.filtered_proteins)
        .join(FILTER_GENE_MODELS.out.filtered_aed)

    RENAME_GENE_MODELS(ch_rename_input)

    // ── Step 6: Get descriptive stats for annotations (per genome) ─────────────
    ch_agat_inputs = ch_assemblies_filtered
        .join(RENAME_GENE_MODELS.out.renamed_gff3, by: 0)
        .map { id, assembly_path, gff_path -> tuple(id, assembly_path, gff_path) }

    AGAT(ch_agat_inputs)

    // ── Step 7: Run BUSCO on protein predictions (per genome) ─────────────
    BUSCO(
        RENAME_GENE_MODELS.out.renamed_proteins,
        DOWNLOAD_BUSCO_DB.out.db)

    // ── Step 8: Run OMArk on protein predictions (per genome) ─────────────
    OMARK(
        RENAME_GENE_MODELS.out.renamed_proteins,
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
        RENAME_GENE_MODELS.out.renamed_proteins)

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

    /*
    ========================================================================================
        Annotation Merging
    ========================================================================================
    */
    // ── Step 17: Merge all functional annotations into single file per genome ─────────────
    ch_func_merged = MERGE_INTERPROSCAN.out.merged_tsv
        .join(MERGE_EGGNOG.out.merged)
        .join(MERGE_DIAMOND.out.merged)
        .join(MERGE_KOFAMSCAN.out.merged)
        .join(MERGE_DEEPSIG.out.merged_gff3)
        .join(RENAME_GENE_MODELS.out.renamed_proteins)
        .join(RENAME_GENE_MODELS.out.renamed_gff3)
        .join(ch_assemblies_filtered)
        .join(OMARK.out.omark_ump)
        .join(RENAME_GENE_MODELS.out.renamed_aed)
        // produces: [assembly_id, ipr_tsv, eggnog_tsv, diamond_tsv, kofam_tsv, deepsig_gff3,
        //            renamed_proteins_fa, renamed_gff3, genome_fasta, omark_ump, renamed_aed]

    MERGE_FUNCTIONAL_ANNOTATIONS(ch_func_merged)

    // ── Step 18: Merge all GFF3s and decorate with functional annotations ───────
    ch_annotate_input = RENAME_GENE_MODELS.out.renamed_gff3
        .join(tRNASCAN.out.trnascan_gff3)
        .join(INFERNAL.out.infernal_gff3)
        .join(TELOMERE_ANNOTATION.out.telomere_gff3)
        .join(CENTROMERE_TO_GFF3.out.centromere_gff3)
        .join(REPEAT_MASK.out.repeat_gff)
        .join(MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv)
    // produces: [assembly_id, braker_gff3, trnascan_gff3, infernal_gff3,
    //            telomere_gff3, centromere_gff3, repeat_gff, func_tsv]

    ANNOTATE_GFF3(ch_annotate_input)

    // ── Step 18b: Structural QC per assembly ─────────────────────────
    ch_structural_qc_in = MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv
        .join(RENAME_GENE_MODELS.out.renamed_gff3)

    STRUCTURAL_QC(ch_structural_qc_in)

    // ── Step 19: Get relevant sequences with annotations ───────
    ch_extract_input = ANNOTATE_GFF3.out.annotated_gff3
        .join(ch_assemblies_filtered, by: 0)
        // produces: [assembly_id, annotated_gff3, filtered_assembly]

    EXTRACT_SEQUENCES(ch_extract_input)

    /*
    ========================================================================================
        Create Report
    ========================================================================================
    */
    ANNOTATION_REPORT(
        FILTER_ASSEMBLY.out.bundle.map { id, folder -> file("${folder}/${id}.assembly_stats.tsv") }.collect(),
        ANNOTATE_GFF3.out.annotated_gff3.map { id, f -> f }.collect(),
        ANNOTATE_GFF3.out.gene_summary.map { id, f -> f }.collect(),
        MERGE_FUNCTIONAL_ANNOTATIONS.out.merged_tsv.map { id, f -> f }.collect(),
        REPEAT_MASK.out.repeat_tbl.map { id, f -> f }.collect(),
        BUSCO.out.busco.map { id, f -> f }.collect(),
        OMARK.out.omark_summary.map { id, f -> f }.collect(),
        OMARK.out.omark_detailed.map { id, f -> f }.collect(),
        AGAT.out.annotation_stats.map { id, f -> f }.collect(),
        tRNASCAN.out.trnascan_tsv.map { id, f -> f }.collect(),
        INFERNAL.out.infernal_tblout.map { id, f -> f }.collect(),
        TELOMERE_ANNOTATION.out.telomere_gff3.map { id, f -> f }.collect(),
        TELOMERE_ANNOTATION.out.telomere_tsv.map { id, f -> f }.collect(),
        CENTROMERE_TO_GFF3.out.centromere_gff3.map { id, f -> f }.collect(),
        STRUCTURAL_QC.out.structural_qc.map { id, f -> f }.collect()
    )
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
