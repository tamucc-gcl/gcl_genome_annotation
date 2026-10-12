#!/usr/bin/env nextflow
// Development entry: accepted selected reads + normalized references -> BRAKER BAMs.
nextflow.enable.dsl = 2
include { RESOLVE_ANNOTATION_INPUTS; ASSERT_ANNOTATION_INPUTS } from './modules/annotationInputs.nf'
include { PLAN_RNA_ALIGNMENT; INDEX_RNA_REFERENCE; ALIGN_RNA_READSET; QC_RNA_ALIGNMENT; RESOLVE_RNA_STRAND; REPORT_RNA_ALIGNMENT; ASSERT_RNA_ALIGNMENT } from './modules/rnaAlignment.nf'

workflow {
    if (!params.sample_sheet) error 'Provide --sample_sheet'
    def encode = { value -> groovy.json.JsonOutput.toJson(value).getBytes('UTF-8').encodeBase64().toString() }
    def request = [launch_dir: launchDir.toString(), sample_sheet: params.sample_sheet.toString(),
        rna_sheet: params.rna_sheet?.toString(), taxdump_dir: params.taxdump_dir?.toString(),
        params: params.collectEntries { key, value -> [key.toString(), value] }]
    def helper = file("${projectDir}/py_scripts/rna_alignment.py", checkIfExists: true)
    def helpers = [helper, file("${projectDir}/py_scripts/reference_repeats.py", checkIfExists: true)]
    def readHelpers = [helper, file("${projectDir}/py_scripts/prepare_rna.py", checkIfExists: true)]
    def prepared = file(params.rna_prepared_status ?: "${params.outdir}/rna/rna_status.json", checkIfExists: true)
    def referenceDir = file(params.rna_reference_dir ?: "${params.outdir}/filtered_assemblies").toString()
    if ((params.rna_min_intron as int) < 1 || (params.rna_max_intron as int) < (params.rna_min_intron as int))
        error 'Invalid RNA intron limits'
    if ((params.rna_strand_min_mapq as int) < 0 || (params.rna_strand_min_mapq as int) > 255)
        error 'RNA strand MAPQ must be between 0 and 255'
    RESOLVE_ANNOTATION_INPUTS(encode(request), file("${projectDir}/py_scripts/resolve_annotation_inputs.py"),
        file("${projectDir}/assets/taxonomy_rules.json"))
    ASSERT_ANNOTATION_INPUTS(RESOLVE_ANNOTATION_INPUTS.out.bundle)
    PLAN_RNA_ALIGNMENT(RESOLVE_ANNOTATION_INPUTS.out.bundle, ASSERT_ANNOTATION_INPUTS.out, prepared,
        encode([reference_dir: referenceDir, strand_policy: [min_informative: params.rna_strand_min_informative,
                dominant_fraction: params.rna_strand_dominant_fraction,
                unstranded_tolerance: params.rna_strand_unstranded_tolerance]]), helper)

    def plans = PLAN_RNA_ALIGNMENT.out.plan.map { path -> new groovy.json.JsonSlurper().parse(path.toFile()) }
    def references = plans.flatMap { plan -> plan.references }
    INDEX_RNA_REFERENCE(references.map { row ->
        tuple(row.reference_id, file(row.reference, checkIfExists: true), file(row.identity, checkIfExists: true))
    }.unique { it[0] }, helpers)
    def jobs = plans.flatMap { plan -> plan.jobs }.map { job -> tuple(job.reference_id, job) }
    def indexed = jobs.combine(INDEX_RNA_REFERENCE.out.index, by: 0).map { refid, job, index ->
        tuple([job.assembly_id, job.readset_id], job, index)
    }
    ALIGN_RNA_READSET(indexed.map { key, job, index ->
        tuple(job.assembly_id, job.readset_id, encode(job.spec), index,
              job.read_paths.collect { file(it, checkIfExists: true) })
    }, readHelpers)
    def aligned = ALIGN_RNA_READSET.out.alignment.map { aid, rid, folder -> tuple([aid, rid], folder) }
    QC_RNA_ALIGNMENT(indexed.join(aligned).map { key, job, index, folder ->
        tuple(job.assembly_id, job.readset_id, encode(job.spec + [min_mapq: params.rna_strand_min_mapq]), folder, index)
    }, helper)
    def qc = QC_RNA_ALIGNMENT.out.qc.map { aid, rid, folder -> tuple([aid, rid], folder) }
    RESOLVE_RNA_STRAND(indexed.join(qc).map { key, job, index, folder ->
        tuple(job.assembly_id, job.readset_id, encode(job.policy), folder)
    }, helper)
    def bams = ALIGN_RNA_READSET.out.alignment.map { aid, rid, folder ->
        [assembly_id: aid, readset_id: rid, bam: "readsets/${aid}/${rid}/aligned.bam",
         bam_index: "readsets/${aid}/${rid}/aligned.bam.csi"]
    }.toList().map { rows -> encode(rows) }
    REPORT_RNA_ALIGNMENT(PLAN_RNA_ALIGNMENT.out.plan,
        QC_RNA_ALIGNMENT.out.qc.map { aid, rid, folder -> file("${folder}/alignment.json", checkIfExists: true) }.toList(),
        RESOLVE_RNA_STRAND.out.strand.map { aid, rid, path -> path }.toList(), bams, helper)
    ASSERT_RNA_ALIGNMENT(REPORT_RNA_ALIGNMENT.out.status)
}
