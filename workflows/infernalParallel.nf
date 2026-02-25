/*
========================================================================================
    INFERNAL PARALLEL — Subworkflow
========================================================================================
    Drop-in replacement for the INFERNAL module that splits the masked genome
    into balanced chunks via greedy bin-packing, runs cmscan in parallel across
    nodes, then merges results.

    Since cmscan reports coordinates relative to original scaffold names and
    whole scaffolds are kept intact (no within-scaffold splitting), the merge
    is a simple concatenation with no coordinate adjustment.

    Inputs  (match original INFERNAL signature):
        - ch_masked_fasta:  tuple val(assembly_id), path(masked_assembly)
        - rfam_database:    path (value channel)

    Outputs (match original INFERNAL signature):
        - infernal_tblout:  tuple val(assembly_id), path("${assembly_id}_rfam.tblout")
        - infernal_cmscan:  tuple val(assembly_id), path("${assembly_id}_rfam.cmscan")
========================================================================================
*/

include { CHUNK_GENOME }    from '../modules/chunk_genome.nf'
include { INFERNAL }        from '../modules/infernal.nf'
include { MERGE_INFERNAL }  from '../modules/infernal.nf'

workflow INFERNAL_PARALLEL {

    take:
    ch_masked_fasta   // tuple val(assembly_id), path(masked_assembly)
    rfam_database     // path (value channel)

    main:

    // ── Step 1: Split genome into balanced chunks ────────────────────────────
    CHUNK_GENOME(ch_masked_fasta)

    // ── Step 2: Flatten chunks into (assembly_id, chunk_id, chunk_file) ──────
    ch_genome_chunks = CHUNK_GENOME.out.chunks
        .flatMap { assembly_id, f_or_list ->
            def files = (f_or_list instanceof Collection) ? f_or_list : [f_or_list]
            files
                .sort { it.name }
                .collect { f ->
                    tuple(assembly_id, f.baseName, f)
                }
        }

    // ── Step 3: Run cmscan per chunk ─────────────────────────────────────────
    INFERNAL(ch_genome_chunks, rfam_database)

    // ── Step 4: Group by assembly and merge ──────────────────────────────────
    ch_tblout = INFERNAL.out.infernal_tblout
        .map { assembly_id, chunk_id, tblout -> tuple(assembly_id, tblout) }
        .groupTuple(size: params.n_infernal_chunks)

    ch_cmscan = INFERNAL.out.infernal_cmscan
        .map { assembly_id, chunk_id, cmscan -> tuple(assembly_id, cmscan) }
        .groupTuple(size: params.n_infernal_chunks)

    ch_merge = ch_tblout
        .join(ch_cmscan, by: 0)

    MERGE_INFERNAL(ch_merge)

    emit:
    infernal_tblout = MERGE_INFERNAL.out.infernal_tblout
    infernal_cmscan = MERGE_INFERNAL.out.infernal_cmscan
}