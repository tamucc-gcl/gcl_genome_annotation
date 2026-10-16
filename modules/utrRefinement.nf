process PREPARE_UTR_INPUTS {
    tag "${aid}"
    label 'aed_python'
    input:
    tuple val(aid), val(enabled), val(readsets), path(gff), path(validation), path(identity), path(rna), path(bams, stageAs: 'bams/readset??.bam'), path(indices, stageAs: 'bams/readset??.bam.csi')
    path helpers
    output:
    tuple val(aid), path('prepared/inputs.json'), emit: status
    script:
    def meta = groovy.json.JsonOutput.toJson(readsets)
    def reads = rna instanceof Collection ? rna : [rna]
    def bamList = bams instanceof Collection ? bams : [bams]
    def indexList = indices instanceof Collection ? indices : [indices]
    def rnaOption = enabled && reads ? "--rna '${reads[0]}'" : ''
    def bamArgs = bamList.collect { "'${it}'" }.join(' ')
    def indexArgs = indexList.collect { "'${it}'" }.join(' ')
    """
    set -euo pipefail
    cat > readsets.json <<'JSON'
    ${meta}
    JSON
    python utr_refinement.py prepare --gff '${gff}' --validation '${validation}' --identity '${identity}' \
      --metadata readsets.json ${rnaOption} --bams ${bamArgs} --indices ${indexArgs} --outdir prepared
    """
}

process RUN_UTR_REFINEMENT {
    tag "${aid}"
    label 'utr_ranking'
    publishDir "${params.outdir}/structural_annotation/utr/${aid}/native", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(enabled), val(readsets), path(gff), path(validation), path(identity), path(rna), path(bams, stageAs: 'bams/readset??.bam'), path(indices, stageAs: 'bams/readset??.bam.csi'), val(approval)
    output:
    tuple val(aid), path('native.gff3'), emit: annotation
    tuple val(aid), path('ranking'), emit: ranking
    script:
    def reads = rna instanceof Collection ? rna : [rna]
    def bamList = bams instanceof Collection ? bams : [bams]
    def fof = readsets.withIndex().collect { row, i -> "${bamList[i]}\t${row.layout == 'paired'}\t${row.strandedness == 'reverse'}" }.join('\n')
    def fofEncoded = (fof ? fof + '\n' : '').getBytes('UTF-8').encodeBase64().toString()
    def nativeCommand = enabled && reads ? """
    awk 'BEGIN{FS=OFS="\\t"} !/^#/ && NF==9 && (\$7=="+" || \$7=="-")' '${reads[0]}' > resolved_rna.gtf
    ingenannot -v 2 -p ${task.cpus} clusterize resolved_rna.gtf ranking/clustered.gtf -g '${gff}' --keep_atts
    ingenannot -v 2 -p ${task.cpus} isoform_ranking ranking/clustered.gtf -f bam.fof \
      --prefix ranking/isoforms --sj_threshold ${params.utr_sj_threshold} \
      --cov_threshold ${params.utr_cov_threshold} --alt_threshold ${params.utr_alt_threshold}
    ingenannot -v 2 utr_refine '${gff}' ranking/isoforms.alternatives.gff native.gff3 --utr_mode rank --onlynew
    """ : "cp '${gff}' native.gff3"
    """
    set -euo pipefail
    mkdir ranking
    python -c "import base64; from pathlib import Path; Path('bam.fof').write_bytes(base64.b64decode('${fofEncoded}'))"
    ${nativeCommand}
    python -c "from importlib.metadata import version; print(version('ingenannot'))" > ranking/version.txt
    """
}

process FINALIZE_UTR_ANNOTATION {
    tag "${aid}"
    label 'aed_python'
    input:
    tuple val(aid), val(skipped), path(original), path(raw_gff), path(identity)
    path helpers
    output:
    tuple val(aid), path("${aid}.utr.gff3"), emit: annotation
    tuple val(aid), path("${aid}.utr.rejections.json"), emit: rejections
    script:
    def command = skipped ? "cp '${original}' '${aid}.utr.gff3'\nprintf '{}\\n' > '${aid}.utr.rejections.json'" :
        "python utr_refinement.py finalize --original '${original}' --native '${raw_gff}' --identity '${identity}' --output '${aid}.utr.gff3'"
    """
    set -euo pipefail
    ${command}
    """
}

process EXTRACT_UTR_SEQUENCES {
    tag "${aid}"
    label 'utr_sequences'
    input:
    tuple val(aid), path(genome), path(original), path(refined)
    output:
    tuple val(aid), path('sequences'), emit: sequences
    script:
    """
    set -euo pipefail
    mkdir sequences
    gffread '${original}' -g '${genome}' -x sequences/baseline.cds.fa -y sequences/baseline.proteins.fa
    gffread '${refined}' -g '${genome}' -w sequences/transcripts.fa -x sequences/cds.fa -y sequences/proteins.fa
    """
}

process VALIDATE_UTR_REFINEMENT {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/utr/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), val(skipped), path(original), path(refined, stageAs: 'candidate.gff3'), path(identity), path(genome), path(inputs), path(sequences), path(rejections)
    path helpers
    output:
    tuple val(aid), path("${aid}.utr.gff3"), emit: annotation
    tuple val(aid), path('utr_status.json'), emit: status
    tuple val(aid), path('utr_changes.tsv'), emit: changes
    tuple val(aid), path('README.md'), emit: report
    tuple val(aid), path('transcripts.fa'), emit: transcripts
    script:
    def skip = skipped ? '--skipped' : ''
    """
    set -euo pipefail
    python utr_refinement.py validate --original '${original}' --refined '${refined}' \
      --identity '${identity}' --genome '${genome}' --inputs '${inputs}' \
      --old-cds '${sequences}/baseline.cds.fa' --new-cds '${sequences}/cds.fa' \
      --old-proteins '${sequences}/baseline.proteins.fa' --new-proteins '${sequences}/proteins.fa' --outdir . ${skip} --rejections '${rejections}'
    cp '${refined}' '${aid}.utr.gff3'
    cp '${sequences}/transcripts.fa' transcripts.fa
    """
}
