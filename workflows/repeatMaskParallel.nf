/*
========================================================================================
    REPEAT MASK PARALLEL — Subworkflow
========================================================================================
*/

include { PREP_RM_BATCHES }     from '../modules/prepRmBatches.nf'
include { REPEAT_MASK_BATCH }   from '../modules/repeatMaskBatch.nf'
include { MERGE_RM_RESULTS }    from '../modules/mergeRmResults.nf'

workflow REPEAT_MASK_PARALLEL {

    take:
    ch_mask_input   // tuple val(assembly_id), path(assembly_path), path(combined_repeats)
    n_batches

    main:

    // At the top of the workflow, create channels for the scripts
    ch_build_script = file("${params.pipeline_root}/py_scripts/build_coord_map.py", checkIfExists: true)
    ch_verify_script = file("${params.pipeline_root}/py_scripts/verify_coord_map.py", checkIfExists: true)
    ch_integrity_script = file("${params.pipeline_root}/py_scripts/reference_repeats.py", checkIfExists: true)

    // ── Split input channel ──────────────────────────────────────────────────
    ch_assemblies = ch_mask_input.map { id, assembly, repeats -> tuple(id, assembly) }
    ch_curated    = ch_mask_input.map { id, assembly, repeats -> tuple(id, repeats) }

    // ── Step 1: Split genome into chunks ─────────────────────────────────────
    PREP_RM_BATCHES(ch_assemblies, n_batches, ch_integrity_script)
    ch_expected = PREP_RM_BATCHES.out.batches.map { aid, files ->
        tuple(aid, files instanceof Collection ? files.size() : 1)
    }

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
    
    ch_batch_full = ch_batch_inputs
        .combine(ch_curated, by: 0)
        .combine(PREP_RM_BATCHES.out.seq_lengths, by: 0)
    
    
    REPEAT_MASK_BATCH(ch_batch_full,
                        ch_build_script,
                        ch_verify_script)

    // ── Step 5: Group by assembly and merge ──────────────────────────────────
    ch_merge = REPEAT_MASK_BATCH.out.batch_results
        .combine(ch_expected, by: 0)
        .map { assembly_id, batch_id, out_file, fasta_file, align_file, count ->
            tuple(groupKey(assembly_id, count), batch_id, out_file, fasta_file, align_file)
        }
        .groupTuple()
        .map { key, batches, outs, fastas, aligns ->
            if (batches.toSet().size() != key.groupSize) error "Duplicate/missing masking batches for ${key.groupTarget}"
            tuple(key.groupTarget, outs, fastas, aligns)
        }
        .join(PREP_RM_BATCHES.out.seq_lengths, by: 0)

    MERGE_RM_RESULTS(ch_merge, ch_integrity_script)

    emit:
    masked_fasta    = MERGE_RM_RESULTS.out.masked_fasta
    repeat_gff      = MERGE_RM_RESULTS.out.repeat_gff
    repeat_tbl      = MERGE_RM_RESULTS.out.repeat_tbl
    repeat_out      = MERGE_RM_RESULTS.out.repeat_out
    repeat_tbl_full = MERGE_RM_RESULTS.out.repeat_tbl_full
    repeat_align    = MERGE_RM_RESULTS.out.repeat_align
    repeat_divsum   = MERGE_RM_RESULTS.out.repeat_divsum
    repeat_landscape = MERGE_RM_RESULTS.out.repeat_landscape
    batch_accounting = PREP_RM_BATCHES.out.accounting
}
