process BUILD_CENTROANNO {
    label 'centromere_build'
    input:
    val revision
    output:
    path('centroanno_tool'), emit: tool
    script:
    if (!(revision ==~ /[0-9a-f]{40}/)) error 'centroAnno requires a pinned Git commit'
    """
    set -euo pipefail
    git init upstream
    git -C upstream fetch --depth 1 https://github.com/junhaiqi/centroAnno.git '${revision}'
    git -C upstream checkout --detach FETCH_HEAD
    git -C upstream submodule update --init --recursive
    make -C upstream -j ${task.cpus} CXX="\${CXX:-g++}"
    mkdir -p centroanno_tool/script
    cp upstream/centroAnno centroanno_tool/
    cp upstream/script/centroFinder.py upstream/script/detectHORs.py centroanno_tool/script/
    git -C upstream rev-parse HEAD > centroanno_tool/revision.txt
    """
    stub:
    """
    mkdir centroanno_tool
    """
}

process CENTROMERE_ANNOTATION {
    tag "${aid}"
    label 'centromere'
    publishDir "${params.outdir}/centromere/${aid}", mode: 'copy', overwrite: true
    input:
    tuple val(aid), path(genome), path(identity), path(tool)
    path helpers
    output:
    tuple val(aid), path("${aid}.centromere.gff3"), emit: centromere_gff3
    path('arrays'), emit: arrays
    path('arrays.log'), emit: runtime_log
    script:
    def threads = Math.min(task.cpus as int, params.centromere_threads as int)
    """
    set -euo pipefail
    python centromere_arrays.py --aid '${aid}' --genome '${genome}' --identity '${identity}' \\
      --tool '${tool}' --threads ${threads} --min-monomer ${params.centromere_min_monomer_length} \\
      --gff '${aid}.centromere.gff3' --outdir arrays 2>&1 | tee arrays.log
    """
    stub:
    """
    mkdir arrays
    printf '##gff-version 3\n' > '${aid}.centromere.gff3'
    touch arrays.log
    """
}
