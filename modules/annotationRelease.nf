process PACKAGE_ANNOTATION_RELEASE {
    tag "${aid}"
    label 'annotation_release'
    publishDir "${params.outdir}/release", mode:'copy', overwrite:true
    input:
    tuple val(aid), path(identified,stageAs:'identified'), path(genome), path(masked), path(identity),
        path(repeats,stageAs:'repeat_source.gff3'), path(repeat_status), path(track_status),
        path(tracks,stageAs:'selected_tracks/*'), val(functional_enabled),
        path(functional_gff,stageAs:'functional/coding.gff3'), path(functional_table,stageAs:'functional/table.tsv'),
        path(functional_genes,stageAs:'functional/genes.tsv'), path(functional_status,stageAs:'functional/status.json')
    path helpers
    output:
    tuple val(aid), path("${aid}"), emit: bundle
    script:
    def trackFiles=tracks instanceof Collection ? tracks : [tracks]
    def request=[aid:aid,identified:identified.toString(),genome:genome.toString(),masked:masked.toString(),identity:identity.toString(),
        repeats:repeats.toString(),repeat_status:repeat_status.toString(),track_status:track_status.toString(),tracks:trackFiles.collect { it.toString() },
        functional_enabled:functional_enabled,functional:functional_enabled ? [gff:functional_gff.toString(),table:functional_table.toString(),
            genes:functional_genes.toString(),status:functional_status.toString()] : [:]]
    def encoded=groovy.json.JsonOutput.toJson(request).getBytes('UTF-8').encodeBase64().toString()
    """
    python package_annotation.py --request '${encoded}' --outdir '${aid}'
    """
}

process ANNOTATION_RELEASE_INDEX {
    label 'annotation_release'
    publishDir "${params.outdir}/release", mode:'copy', overwrite:true
    input:
    val assemblies
    path bundles
    output:
    path 'README.md'
    path 'release_index.json'
    script:
    def encoded=groovy.json.JsonOutput.toJson(assemblies).getBytes('UTF-8').encodeBase64().toString()
    """
    python - <<'PY'
    import base64,json
    from pathlib import Path
    assemblies=json.loads(base64.b64decode('${encoded}'))
    rows=[json.loads((Path(a)/'release_status.json').read_text()) for a in assemblies]
    assert all(r['status']=='ready' and r['assembly_id']==a for a,r in zip(assemblies,rows))
    Path('release_index.json').write_text(json.dumps(dict(status='ready',assemblies=rows),indent=2)+'\\n')
    lines=['# Final annotation releases','','| Assembly | Coding genes | Coding transcripts | Deliverables |','|---|---:|---:|---|']
    lines += ['| '+a+' | '+str(r['genes'])+' | '+str(r['coding_transcripts'])+' | [Release]('+a+'/README.md) |' for a,r in zip(assemblies,rows)]
    lines += ['','Each release manifest lists the current deliverables and checksums. Raw results and detailed stage reports remain in their original directories. Earlier published files may persist after settings change; the current manifest is authoritative.','']
    Path('README.md').write_text('\\n'.join(lines))
    PY
    """
}
