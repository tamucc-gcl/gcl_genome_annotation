/*
========================================================================================
    CURATE DFAM
========================================================================================
    Extract taxon-specific repeat library from Dfam database
    Runs ONCE per pipeline execution (not per genome)
    
    Inputs:
        - dfam_dir:   Path to Dfam database directory
        - taxa_name:  Taxonomic group name or NCBI taxid (e.g. 'Actinopterygii')
    
    Outputs:
        - dfam_lib: dfam_subset.lib (FASTA-format repeat library)
========================================================================================
*/

process CURATE_DFAM {
    tag "${taxa_name}"
    label 'dfam'
    storeDir "${params.dfam_database}/subsets/dfam${params.dfam_version}/tetools-${params.tetools_version}"
    input:
    tuple val(source_key), val(taxa_name), path(dfam_dir)
    output:
    tuple val(source_key), path("${source_key}"), emit: bundle
    script:
    """
    set -euo pipefail
    mkdir -p ${source_key}
    famdb.py -i '${dfam_dir}' info > ${source_key}/famdb_info.txt
    famdb.py -i '${dfam_dir}' families --format fasta_name --include-class-in-name \
        --ancestors --descendants '${taxa_name}' > ${source_key}/dfam_subset.lib
    python - <<'PY'
    import hashlib, json, re
    from pathlib import Path
    folder = Path('${source_key}')
    Path('${params.dfam_database}/subsets/dfam${params.dfam_version}/tetools-${params.tetools_version}').mkdir(parents=True, exist_ok=True)
    info = (folder / 'famdb_info.txt').read_text()
    expected = '${params.dfam_version}'
    version = expected[0] + '.' + expected[1:]
    if not re.search(r'Version\\s*:\\s*' + re.escape(version) + r'(?:\\s|\\Z)', info):
        raise SystemExit('FamDB release does not match the configured Dfam release: ' + version)
    library = folder / 'dfam_subset.lib'
    if not library.stat().st_size:
        raise SystemExit('Dfam extraction produced an empty library for ${taxa_name}')
    files = sorted(Path('${dfam_dir}').glob('dfam' + expected + '_full.*.h5'))
    partitions = {int(path.name.split('.')[-2]) for path in files}
    if partitions != set(range(int('${params.dfam_partition_count}'))):
        raise SystemExit('Incomplete/mixed FamDB partitions for the configured full release')
    (folder / 'dfam_subset.source.json').write_text(json.dumps({
        'release': version, 'taxa_name': '${taxa_name}', 'tetools': '${params.tetools_version}',
        'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(),
        'partitions': [{'name': path.name, 'bytes': path.stat().st_size} for path in files],
        'source_identity_scope': 'Declared release, FamDB info and partition sizes; subset is checksummed'
    }, indent=2) + '\\n')
    PY
    """
}

/*
========================================================================================
    DOWNLOAD DFam DATABASE MODULE
========================================================================================
*/

process DOWNLOAD_DFAM_DB {
    tag "${dfam_version}"
    label 'download_dfam_db'
    storeDir "${params.dfam_database}"

    input:
    val dfam_version

    output:
    path("dfam_db"), emit: dfam_db_dir

    script:
    """
    set -euo pipefail
    mkdir -p dfam_db
    cd dfam_db

    wget -q -O - ${params.dfam_release_url}/ \
        | grep -oP "dfam${dfam_version}_full\\..*?\\.h5\\.gz" \
        | sort -u \
        | awk '{print "${params.dfam_release_url}/" \$0}' \
        > files.txt

    cpus=${task.cpus}

    # Conservative split
    parallel_files=\$(( cpus / 8 ))
    if [ \${parallel_files} -lt 1 ]; then
        parallel_files=1
    fi
    connections_per_file=8

    aria2c \
        -j \${parallel_files} \
        -x \${connections_per_file} \
        -s \${connections_per_file} \
        -i files.txt

    pigz -p ${task.cpus} -d "dfam${dfam_version}_full."*.h5.gz
    """
}
