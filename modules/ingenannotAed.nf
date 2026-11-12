process PREPARE_INGENANNOT_EVIDENCE {
    tag "${aid}"
    label 'aed_index'
    input:
    tuple val(aid), path(gff), path(validation), path(rna, stageAs: 'rna/*'), path(protein, stageAs: 'protein/*'), path(provenance, stageAs: 'provenance/*')
    path helpers
    output:
    tuple val(aid), path('indexed'), emit: evidence
    script:
    def r = rna instanceof Collection ? rna : [rna]
    def p = protein instanceof Collection ? protein : [protein]
    def meta = provenance instanceof Collection ? provenance : [provenance]
    def options = (r ? "--rna '${r[0]}'" : '') + ' ' + (p ? "--protein '${p[0]}' --protein-provenance '${meta[0]}'" : '')
    def rpath = r ? "'${r[0]}'" : '/dev/null'
    """
    set -euo pipefail
    mkdir indexed
    python ingenannot_inputs.py --gff '${gff}' --validation '${validation}' ${options} --split-protein-pools --out indexed/inputs.json
    LC_ALL=C sort -k1,1 -k4,4n -k5,5n ${rpath} | bgzip -c > indexed/rna.gtf.gz
    tabix -p gff indexed/rna.gtf.gz
    for pool in reference orthodb; do
      LC_ALL=C sort -k1,1 -k4,4n -k5,5n indexed/\$pool.gff | bgzip -c > indexed/\$pool.gff.gz
      tabix -p gff indexed/\$pool.gff.gz
      rm indexed/\$pool.gff
    done
    """
}

process RUN_INGENANNOT_AED {
    tag "${aid}"
    label 'ingenannot'
    publishDir "${params.outdir}/structural_annotation/aed/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(gff), path(indexed)
    output:
    tuple val(aid), path('native'), emit: annotated
    script:
    def cds = params.ingenannot_rna_cds_only ? '--aed_tr_cds_only' : ''
    """
    set -euo pipefail
    mkdir native
    for pool in reference orthodb; do
      if [ "\$pool" = orthodb ] && ! python -c 'import json,sys; sys.exit(not json.load(open(sys.argv[1]))["evidence"]["protein"]["available"])' '${indexed}/orthodb.json'; then
        continue
      fi
      ingenannot -v 2 -p ${task.cpus} aed '${gff}' native/\$pool.aed.gff3 braker3 \\
        '${indexed}/rna.gtf.gz' "${indexed}/\$pool.gff.gz" \\
        --evtr_source stringtie --evpr_source miniprot --evtrstranded --evprstranded \\
        --penalty_overflow ${params.ingenannot_rna_penalty} ${cds}
    done
    python -c "from importlib.metadata import version; print(version('ingenannot'))" > native/ingenannot_version.txt
    """
}

