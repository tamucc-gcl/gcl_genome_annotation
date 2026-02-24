/*
========================================================================================
    REPEAT MASK PARALLEL — Subworkflow
========================================================================================
    Drop-in replacement for the REPEAT_MASK module that splits the genome into
    chunks, runs RepeatMasker in parallel, then merges and adjusts results.

    Interface matches the original REPEAT_MASK module:

    Inputs:
        - ch_mask_input:  tuple val(assembly_id), path(assembly_path), path(curated_repeats)
        - dfam_lib:       path (Dfam taxon-specific repeat library)

    Outputs:
        - masked_fasta:    tuple val(assembly_id), path("${assembly_id}.masked.fasta")
        - repeat_gff:      tuple val(assembly_id), path("${assembly_id}.gff")
        - repeat_tbl:      tuple val(assembly_id), path("${assembly_id}.tbl")
        - repeat_out:      tuple val(assembly_id), path("${assembly_id}.out")
        - repeat_tbl_full: tuple val(assembly_id), path("${assembly_id}.buildSummary.tbl")
========================================================================================
*/

include { PREP_RM_BATCHES }     from '../modules/prepRmBatches.nf'
include { COMBINE_REPEAT_LIBS } from '../modules/combineRepeatLibs.nf'
include { REPEAT_MASK_BATCH }   from '../modules/repeatMaskBatch.nf'
include { MERGE_RM_RESULTS }    from '../modules/mergeRmResults.nf'

workflow REPEAT_MASK_PARALLEL {

    take:
    ch_mask_input   // tuple val(assembly_id), path(assembly_path), path(curated_repeats)
    dfam_lib        // path (value channel)

    main:

    // ── Split input channel ──────────────────────────────────────────────────
    ch_assemblies = ch_mask_input.map { id, assembly, repeats -> tuple(id, assembly) }
    ch_curated    = ch_mask_input.map { id, assembly, repeats -> tuple(id, repeats) }

    // ── Step 1: Split genome into chunks ─────────────────────────────────────
    PREP_RM_BATCHES(ch_assemblies, params.n_rm_chunks)

    // ── Step 2: Combine repeat libraries (once per assembly) ─────────────────
    COMBINE_REPEAT_LIBS(ch_curated, dfam_lib)

    // ── Step 3: Flatten batches and pair with BEDs ───────────────────────────
    ch_batches = PREP_RM_BATCHES.out.batches
        .flatMap { assembly_id, fa_list ->
            def files = (fa_list instanceof Collection) ? fa_list : [fa_list]
            files.sort { it.name }.collect { f ->
                tuple(assembly_id, f.baseName, f)
            }
        }

    ch_beds = PREP_RM_BATCHES.out.beds
        .flatMap { assembly_id, bed_list ->
            def files = (bed_list instanceof Collection) ? bed_list : [bed_list]
            files.sort { it.name }.collect { f ->
                tuple(assembly_id, f.baseName, f)
            }
        }

    ch_batch_inputs = ch_batches.join(ch_beds, by: [0, 1])

    // ── Step 4: Run RepeatMasker per batch + adjust coordinates ──────────────
    REPEAT_MASK_BATCH(
        ch_batch_inputs,
        COMBINE_REPEAT_LIBS.out.combined_lib.map { id, lib -> lib }.collect(),
        PREP_RM_BATCHES.out.seq_lengths.map { id, f -> f }.collect()
    )

    // ── Step 5: Group by assembly and merge ──────────────────────────────────
    ch_merge = REPEAT_MASK_BATCH.out.batch_results
        .map { assembly_id, batch_id, out_file, fasta_file ->
            tuple(assembly_id, out_file, fasta_file)
        }
        .groupTuple(size: params.n_rm_chunks)

    MERGE_RM_RESULTS(ch_merge)

    emit:
    masked_fasta    = MERGE_RM_RESULTS.out.masked_fasta
    repeat_gff      = MERGE_RM_RESULTS.out.repeat_gff
    repeat_tbl      = MERGE_RM_RESULTS.out.repeat_tbl
    repeat_out      = MERGE_RM_RESULTS.out.repeat_out
    repeat_tbl_full = MERGE_RM_RESULTS.out.repeat_tbl_full
}