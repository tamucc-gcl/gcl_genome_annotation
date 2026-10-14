// Existing miniprot alignment command separated from the legacy scoring process.
process ALIGN_AED_PROTEINS {
    tag "${aid}"
    label 'computeAed'
    publishDir "${params.outdir}/structural_annotation/evidence_aed/${aid}/protein_alignment", mode: params.publish_dir_mode
    input:
    tuple val(aid), path(genome), path(proteins), path(identity), path(provenance)
    path helpers
    output:
    tuple val(aid), path("${aid}.miniprot.gff"), emit: alignment
    tuple val(aid), path('miniprot_version.txt'), emit: version
    script:
    """
    set -euo pipefail
    python - <<'PY'
    import hashlib, json
    from pathlib import Path
    from reference_repeats import fasta
    from prepare_protein_evidence import checksum
    identity = json.loads(Path('${identity}').read_text())
    contigs = [{'id': header.split()[0], 'length': len(sequence),
               'uppercase_sha256': hashlib.sha256(sequence.upper().encode('ascii')).hexdigest()}
              for header, sequence in fasta('${genome}')]
    if contigs != identity['contigs']:
        raise SystemExit('Protein alignment genome differs from accepted reference')
    source = json.loads(Path('${provenance}').read_text())
    if source['status'] != 'ready' or checksum('${proteins}') != source['proteins_sha256']:
        raise SystemExit('Protein evidence differs from accepted source manifest')
    PY
    miniprot -t ${task.cpus} --gff '${genome}' '${proteins}' > '${aid}.miniprot.gff'
    miniprot --version > miniprot_version.txt
    """
}

process ASSEMBLE_AED_RNA {
    tag "${aid}/${rid}"
    label 'aed_stringtie'
    input:
    tuple val(aid), val(rid), val(strand), path(bam), path(bai)
    output:
    tuple val(aid), path("${rid}.gtf"), emit: transcripts
    script:
    def stranded = strand == 'reverse' ? '--rf' : strand == 'forward' ? '--fr' : ''
    """
    set -euo pipefail
    stringtie '${bam}' -p ${task.cpus} ${stranded} -o '${rid}.gtf'
    """
}

process MERGE_AED_RNA {
    tag "${aid}"
    label 'aed_stringtie'
    publishDir "${params.outdir}/structural_annotation/evidence_aed/${aid}/rna_assembly", mode: params.publish_dir_mode
    input:
    tuple val(aid), path(transcripts, stageAs: 'transcripts/*')
    output:
    tuple val(aid), path('transcripts_merged.gtf'), emit: transcripts
    script:
    def files = transcripts instanceof Collection ? transcripts : [transcripts]
    def arguments = files.collect { "'${it}'" }.join(' ')
    """
    set -euo pipefail
    stringtie --merge -o transcripts_merged.gtf ${arguments}
    """
}

process NORMALIZE_AED_EVIDENCE {
    tag "${aid}"
    label 'aed_python'
    input:
    tuple val(aid), path(identity), path(rna, stageAs: 'rna/*'), path(protein, stageAs: 'protein/*')
    path helpers
    output:
    tuple val(aid), path('evidence'), emit: evidence
    script:
    def asList = { v -> v instanceof Collection ? v : [v] }
    def r = asList(rna)
    def p = asList(protein)
    def options = (r ? "--rna '${r[0]}'" : '') + ' ' + (p ? "--protein '${p[0]}'" : '')
    """
    python evidence_aed.py normalize --identity '${identity}' ${options} --outdir evidence
    """
}

process SCORE_EVIDENCE_AED {
    tag "${aid}"
    label 'aed_python'
    publishDir "${params.outdir}/structural_annotation/evidence_aed/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(models), path(identity), path(evidence), path(validation)
    path helpers
    output:
    tuple val(aid), path('aed_scores.tsv'), emit: scores
    tuple val(aid), path('aed_status.json'), emit: status
    tuple val(aid), path('README.md'), emit: report
    script:
    """
    python evidence_aed_selection.py --models '${models}' --identity '${identity}' \
        --normalized '${evidence}' --validation '${validation}' --selection '${params.aed_selection}' --outdir .
    """
}
