/*
========================================================================================
    REPEAT MASK PARALLEL — Subworkflow
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
    //skip real repeat library for testing
    /*
    ch_batch_full = ch_batch_inputs
        .combine(COMBINE_REPEAT_LIBS.out.combined_lib, by: 0)
        .combine(PREP_RM_BATCHES.out.seq_lengths, by: 0)
    */
    ch_batch_full = ch_batch_inputs
        .combine(ch_curated, by: 0)
        .combine(PREP_RM_BATCHES.out.seq_lengths, by: 0)

    REPEAT_MASK_BATCH(ch_batch_full)

    // ── Step 5: Group by assembly and merge ──────────────────────────────────
    ch_merge = REPEAT_MASK_BATCH.out.batch_results
        .map { assembly_id, batch_id, out_file, fasta_file, align_file ->
            tuple(assembly_id, out_file, fasta_file, align_file)
        }
        .groupTuple(size: params.n_rm_chunks)

    MERGE_RM_RESULTS(ch_merge)

    emit:
    masked_fasta    = MERGE_RM_RESULTS.out.masked_fasta
    repeat_gff      = MERGE_RM_RESULTS.out.repeat_gff
    repeat_tbl      = MERGE_RM_RESULTS.out.repeat_tbl
    repeat_out      = MERGE_RM_RESULTS.out.repeat_out
    repeat_tbl_full = MERGE_RM_RESULTS.out.repeat_tbl_full
    repeat_align    = MERGE_RM_RESULTS.out.repeat_align
    repeat_divsum   = MERGE_RM_RESULTS.out.repeat_divsum
    repeat_landscape = MERGE_RM_RESULTS.out.repeat_landscape
}