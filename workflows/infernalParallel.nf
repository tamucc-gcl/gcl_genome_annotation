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

    Chunk count is determined dynamically from CHUNK_GENOME output rather than
    hardcoded from params.n_infernal_chunks. This handles assemblies with fewer
    scaffolds than requested chunks (greedy bin-packing assigns whole scaffolds,
    so empty bins produce no output file). groupKey() enables per-assembly
    early release so each assembly's merge starts as soon as its chunks finish.

    Inputs  (match original INFERNAL signature):
        - ch_masked_fasta:  tuple val(assembly_id), path(masked_assembly)
        - rfam_database:    path (value channel)

    Outputs (match original INFERNAL signature):
        - infernal_tblout:  tuple val(assembly_id), path("${assembly_id}_rfam.tblout")
        - infernal_cmscan:  tuple val(assembly_id), path("${assembly_id}_rfam.cmscan")
        - infernal_gff3:    tuple val(assembly_id), path("${assembly_id}_rfam.gff3")
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

    // ── Step 1: Split genome into balanced chunks ────────────────────────
    CHUNK_GENOME(ch_masked_fasta)

    // ── Step 2: Count actual chunks per assembly ─────────────────────────
    //   May be fewer than params.n_infernal_chunks if the assembly has
    //   fewer scaffolds than requested bins (greedy bin-packing assigns
    //   whole scaffolds, so empty bins produce no output file).
    ch_chunk_count = CHUNK_GENOME.out.chunks
        .map { assembly_id, files ->
            def fileList = (files instanceof Collection) ? files : [files]
            tuple(assembly_id, fileList.size())
        }

    // ── Step 3: Flatten chunks into (assembly_id, chunk_id, chunk_file) ──
    ch_genome_chunks = CHUNK_GENOME.out.chunks
        .flatMap { assembly_id, f_or_list ->
            def files = (f_or_list instanceof Collection) ? f_or_list : [f_or_list]
            files
                .sort { it.name }
                .collect { f ->
                    tuple(assembly_id, f.baseName, f)
                }
        }

    // ── Step 4: Run cmscan per chunk ─────────────────────────────────────
    INFERNAL(ch_genome_chunks, rfam_database)

    // ── Step 5: Join per-chunk results, group with dynamic size, merge ───
    //   First join all three output types per (assembly_id, chunk_id),
    //   then combine with the actual chunk count so groupKey() enables
    //   per-assembly early release in groupTuple().
    ch_per_chunk = INFERNAL.out.infernal_tblout
        .join(INFERNAL.out.infernal_cmscan, by: [0, 1])
        .join(INFERNAL.out.infernal_gff3, by: [0, 1])
        // tuple(assembly_id, chunk_id, tblout, cmscan, gff3)
        .map { assembly_id, chunk_id, tblout, cmscan, gff3 ->
            tuple(assembly_id, tblout, cmscan, gff3)
        }

    ch_merge = ch_per_chunk
        .combine(ch_chunk_count, by: 0)
        // tuple(assembly_id, tblout, cmscan, gff3, n_chunks)
        .map { assembly_id, tblout, cmscan, gff3, n ->
            tuple(groupKey(assembly_id, n), tblout, cmscan, gff3)
        }
        .groupTuple()
        // tuple(assembly_id, [tblout_files], [cmscan_files], [gff3_files])

    MERGE_INFERNAL(ch_merge)

    emit:
    infernal_tblout = MERGE_INFERNAL.out.infernal_tblout
    infernal_cmscan = MERGE_INFERNAL.out.infernal_cmscan
    infernal_gff3   = MERGE_INFERNAL.out.infernal_gff3
}