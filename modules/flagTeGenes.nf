// Adapted from FLAG_TE_GENES: retain SeqKit/DIAMOND/BEDTools scientific tools,
// separate cached search from threshold-dependent flagging, never remove models.
process BUILD_TE_REPEAT_DB {
    tag "${aid}"
    label 'flag_te_genes'
    publishDir "${params.outdir}/structural_annotation/assessment/${aid}/te_database", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(library)
    path helpers
    output:
    tuple val(aid), path('database'), emit: database
    script:
    """
    set -euo pipefail
    mkdir database
    if [ -s '${library}' ]; then
      seqkit translate -f 6 -F --clean '${library}' | seqkit seq -m 30 > database/repeats.faa
    else
      touch database/repeats.faa
    fi
    if [ -s database/repeats.faa ]; then
      diamond makedb --in database/repeats.faa --db database/repeats.dmnd --threads ${task.cpus}
    fi
    seqkit version > database/seqkit_version.txt
    diamond version > database/diamond_version.txt
    python te_search.py database --library '${library}' --translated database/repeats.faa --folder database
    """
}

process SEARCH_TE_PROTEINS {
    tag "${aid}"
    label 'flag_te_genes'
    publishDir "${params.outdir}/structural_annotation/assessment/${aid}/te_search", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(proteins), path(database)
    path helpers
    output:
    tuple val(aid), path('search'), emit: search
    script:
    """
    set -euo pipefail
    python te_search.py search --proteins '${proteins}' --db '${database}' --threads ${task.cpus} --evalue ${params.te_flag_evalue} --outdir search
    diamond version > search/diamond_version.txt
    """
}

process FLAG_TE_GENES {
    tag "${aid}"
    label 'flag_te_genes'
    publishDir "${params.outdir}/structural_annotation/assessment/${aid}/te", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(gff), path(metrics), path(regions), path(identity), path(repeats), path(search)
    path helpers
    output:
    tuple val(aid), path("${aid}.te_flagged.gff3"), emit: flagged_gff3
    tuple val(aid), path('te_gene_report.tsv'), emit: te_gene_report
    tuple val(aid), path('te_status.json'), emit: status
    script:
    """
    set -euo pipefail
    python te_evidence.py repeat-bed --source '${repeats}' --identity '${identity}' --out repeats.bed
    if [ -s repeats.bed ]; then
      bedtools sort -i repeats.bed | bedtools merge -i - > repeats.union.bed
    else
      touch repeats.union.bed
    fi
    bedtools coverage -a '${regions}' -b repeats.union.bed > coverage.tsv
    python te_evidence.py flag --assembly-id '${aid}' --gff '${gff}' --metrics '${metrics}' --coverage coverage.tsv --search-dir '${search}' --overlap ${params.te_flag_overlap_frac} --outdir .
    """
}
