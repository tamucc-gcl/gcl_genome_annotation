/*
========================================================================================
    FILTER REPEAT LIBRARY
========================================================================================
    Remove gene-derived sequences from de novo repeat library using a
    two-database BLASTX strategy:

      1. BLASTX against RepeatPeps.lib (curated TE proteins from RepeatMasker)
         → any hit = confirmed TE protein → KEEP
      2. BLASTX against a general protein database (e.g. SwissProt)
         → hits that did NOT match RepeatPeps = likely gene-derived → REMOVE

    Repeats with no protein hits at all are retained (non-coding repeats,
    satellites, SINEs, etc.).

    This approach is taxa-independent: RepeatPeps.lib covers all known TE
    superfamilies across eukaryotes, so no keyword whitelist is needed.

    Inputs:
        - assembly_id:       Sample identifier
        - repeat_lib:        De novo curated/classified repeat library FASTA
        - diamond_db:        General protein DIAMOND database (e.g. SwissProt .dmnd)

    Outputs:
        - filtered_lib:      ${assembly_id}_repeats_filtered.fa  (cleaned library)
        - removed_lib:       ${assembly_id}_repeats_removed.fa   (gene-derived seqs)
        - filter_report:     ${assembly_id}_repeat_filter_report.tsv

    Tools: DIAMOND (blastx, makedb), Python 3 (stdlib only), seqkit
    Container: dfam/tetools (for RepeatPeps.lib)
========================================================================================
*/

/*
========================================================================================
    EXTRACT REPEAT PEPS
========================================================================================
    Copy RepeatPeps.lib from the TEtools container to a pipeline-accessible
    path. Runs once per pipeline execution.

    RepeatPeps.lib is a curated collection of TE-encoded protein sequences
    covering all major transposon superfamilies across eukaryotes. It ships
    with RepeatMasker inside the Dfam TEtools container.

    Inputs: none (uses container path)
    Outputs: repeat_peps — RepeatPeps.lib FASTA
========================================================================================
*/

process EXTRACT_REPEAT_PEPS {
    tag "RepeatPeps"
    label 'dfam'
    //storeDir "${params.outdir}/repeats/databases"
    storeDir "${params.dfam_database}"

    output:
    path("RepeatPeps.lib"), emit: repeat_peps

    script:
    """
    cp /opt/RepeatMasker/Libraries/RepeatPeps.lib RepeatPeps.lib
    """

    stub:
    """
    touch RepeatPeps.lib
    """
}


// Share the derived TE-protein index across assemblies/projects, keyed by input content.
process INDEX_REPEAT_PEPS {
    tag "${source_key}"
    label 'repeat_peps_index'
    storeDir "${params.dfam_database}/indices/diamond-${params.repeat_diamond_version}"
    input:
    tuple val(source_key), path(repeat_peps)
    output:
    tuple val(source_key), path("${source_key}.dmnd"), path("${source_key}.source.json"), emit: index
    script:
    """
    set -euo pipefail
    diamond makedb --in '${repeat_peps}' --db '${source_key}' --threads ${task.cpus}
    python - <<'PY'
    import hashlib, json
    from pathlib import Path
    Path('${params.dfam_database}/indices/diamond-${params.repeat_diamond_version}').mkdir(parents=True, exist_ok=True)
    Path('${source_key}.source.json').write_text(json.dumps({
        'input_sha256': hashlib.sha256(Path('${repeat_peps}').read_bytes()).hexdigest(),
        'diamond_version': '${params.repeat_diamond_version}',
        'index_sha256': hashlib.sha256(Path('${source_key}.dmnd').read_bytes()).hexdigest()
    }, indent=2) + '\\n')
    PY
    """
}

process FILTER_REPEAT_LIBRARY {
    tag "${assembly_id}"
    label 'filter_repeat_lib'
    publishDir "${params.outdir}/repeats/denovo", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(repeat_lib)
    path diamond_db
    path repeat_peps
    path repeat_peps_db

    output:
    tuple val(assembly_id), path("${assembly_id}_repeats_filtered.fa"),
        emit: filtered_lib
    tuple val(assembly_id), path("${assembly_id}_repeats_removed.fa"),
        emit: removed_lib
    tuple val(assembly_id), path("${assembly_id}_repeat_filter_report.tsv"),
        emit: filter_report
    tuple val(assembly_id), path("${assembly_id}.repeat_vs_tepeps.tsv"), path("${assembly_id}.repeat_vs_prot.tsv"), emit: hits
    tuple val(assembly_id), path("${assembly_id}.repeat_filter_sources.json"), emit: provenance

    script:
    def max_evalue = params.repeat_filter_evalue != null ? params.repeat_filter_evalue : 1e-10
    def min_pident = params.repeat_filter_pident != null ? params.repeat_filter_pident : 25
    """
    set -euo pipefail
    if [ ! -s '${repeat_lib}' ]; then
        touch repeat_vs_tepeps.tsv repeat_vs_prot.tsv
    else

    # ── 2. BLASTX repeats against RepeatPeps (TE protein database) ────────────
    diamond blastx \
        --db ${repeat_peps_db} \
        --query ${repeat_lib} \
        --out repeat_vs_tepeps.tsv \
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \
        --evalue ${max_evalue} \
        --max-target-seqs 1 \
        --threads ${task.cpus}

    # ── 3. BLASTX repeats against general protein DB (e.g. SwissProt) ─────────
    diamond blastx \
        --db ${diamond_db} \
        --query ${repeat_lib} \
        --out repeat_vs_prot.tsv \
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \
        --evalue ${max_evalue} \
        --max-target-seqs 1 \
        --threads ${task.cpus}
    fi

    # ── 4. Classify using two-database logic ──────────────────────────────────
    python3 - <<'PYEOF'
import sys, hashlib, json
from pathlib import Path

MIN_PIDENT = ${min_pident}

# ── Parse RepeatPeps hits (TE proteins) ───────────────────────────────────
te_hits = {}  # qseqid → best hit
with open("repeat_vs_tepeps.tsv") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue
        qid = cols[0]
        bitscore = float(cols[8])
        if qid not in te_hits or bitscore > te_hits[qid]["bitscore"]:
            te_hits[qid] = {
                "sseqid": cols[1], "stitle": cols[2],
                "pident": float(cols[3]), "evalue": float(cols[7]),
                "bitscore": bitscore,
            }

# ── Parse general protein DB hits ─────────────────────────────────────────
prot_hits = {}
with open("repeat_vs_prot.tsv") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue
        qid = cols[0]
        bitscore = float(cols[8])
        if qid not in prot_hits or bitscore > prot_hits[qid]["bitscore"]:
            prot_hits[qid] = {
                "sseqid": cols[1], "stitle": cols[2],
                "pident": float(cols[3]), "evalue": float(cols[7]),
                "bitscore": bitscore,
            }

# ── Classification logic ─────────────────────────────────────────────────
#   1. Hit RepeatPeps                    → KEEP (confirmed TE protein)
#   2. Hit general DB only, pident >= N  → REMOVE (gene-derived)
#   3. Hit general DB only, pident < N   → KEEP (weak/ambiguous hit)
#   4. No protein hits                   → KEEP (non-coding repeat)

remove_ids = set()
report_rows = []

# All repeat IDs that had any hit
all_query_ids = set(te_hits.keys()) | set(prot_hits.keys())
library_ids = {line[1:].split()[0] for line in Path('${repeat_lib}').read_text().splitlines() if line.startswith('>')}
if not all_query_ids <= library_ids:
    raise ValueError('DIAMOND hits reference an unknown repeat-library ID')
all_query_ids |= library_ids

for qid in sorted(all_query_ids):
    has_te = qid in te_hits
    has_prot = qid in prot_hits

    te = te_hits.get(qid, {})
    pr = prot_hits.get(qid, {})

    if has_te:
        # RepeatPeps hit → confirmed TE, keep regardless of general DB hit
        verdict = "KEEP_te_protein"
        report_rows.append((
            qid,
            te.get("sseqid", "-"), f"{te.get('pident', 0):.1f}",
            f"{te.get('evalue', 0):.2e}", f"{te.get('bitscore', 0):.1f}",
            pr.get("sseqid", "-"), f"{pr.get('pident', 0):.1f}" if has_prot else "-",
            f"{pr.get('evalue', 0):.2e}" if has_prot else "-",
            pr.get("stitle", "-") if has_prot else "-",
            verdict,
        ))
    elif has_prot:
        # General DB hit only
        if pr["pident"] >= MIN_PIDENT:
            verdict = "REMOVE_gene_derived"
            remove_ids.add(qid)
        else:
            verdict = "KEEP_below_threshold"

        report_rows.append((
            qid,
            "-", "-", "-", "-",
            pr["sseqid"], f"{pr['pident']:.1f}",
            f"{pr['evalue']:.2e}",
            pr["stitle"],
            verdict,
        ))

    else:
        report_rows.append((qid, '-', '-', '-', '-', '-', '-', '-', '-', 'KEEP_no_protein_hit'))

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

Path('${assembly_id}.repeat_filter_sources.json').write_text(json.dumps({
    'repeat_library_sha256': sha('${repeat_lib}'), 'general_database_sha256': sha('${diamond_db}'),
    'repeat_peps_sha256': sha('${repeat_peps}'), 'minimum_identity': MIN_PIDENT,
    'maximum_evalue': ${max_evalue}}, indent=2) + '\\n')

# ── Write report ──────────────────────────────────────────────────────────
with open("${assembly_id}_repeat_filter_report.tsv", "w") as out:
    out.write("repeat_id\\t"
              "te_hit_id\\tte_pident\\tte_evalue\\tte_bitscore\\t"
              "prot_hit_id\\tprot_pident\\tprot_evalue\\tprot_description\\t"
              "verdict\\n")
    for row in report_rows:
        out.write("\\t".join(str(x) for x in row) + "\\n")

# ── Write ID list for seqkit ──────────────────────────────────────────────
with open("remove_ids.txt", "w") as out:
    for rid in sorted(remove_ids):
        out.write(rid + "\\n")

# ── Summary ───────────────────────────────────────────────────────────────
n_te = sum(1 for r in report_rows if r[-1] == "KEEP_te_protein")
n_below = sum(1 for r in report_rows if r[-1] == "KEEP_below_threshold")
n_remove = len(remove_ids)
print(f"Repeat library filter (two-database strategy):", file=sys.stderr)
print(f"  RepeatPeps hit (confirmed TE, kept): {n_te}", file=sys.stderr)
print(f"  General DB only, below threshold (kept): {n_below}", file=sys.stderr)
print(f"  General DB only, gene-derived (removed): {n_remove}", file=sys.stderr)

PYEOF

    # ── 5. Split library into kept and removed ────────────────────────────────
    if [ -s remove_ids.txt ]; then
        seqkit grep -v -f remove_ids.txt ${repeat_lib} > ${assembly_id}_repeats_filtered.fa
        seqkit grep -f remove_ids.txt ${repeat_lib} > ${assembly_id}_repeats_removed.fa
    else
        cp ${repeat_lib} ${assembly_id}_repeats_filtered.fa
        touch ${assembly_id}_repeats_removed.fa
    fi

    # ── 6. Summary ────────────────────────────────────────────────────────────
    mv repeat_vs_tepeps.tsv ${assembly_id}.repeat_vs_tepeps.tsv
    mv repeat_vs_prot.tsv ${assembly_id}.repeat_vs_prot.tsv
    n_before=\$(grep -c '^>' ${repeat_lib} || true)
    n_after=\$(grep -c '^>' ${assembly_id}_repeats_filtered.fa || true)
    n_removed=\$((n_before - n_after))
    echo "Repeat library: \${n_before} → \${n_after} sequences (\${n_removed} gene-derived removed)" >&2
    """

    stub:
    """
    touch ${assembly_id}_repeats_filtered.fa
    touch ${assembly_id}_repeats_removed.fa
    touch ${assembly_id}_repeat_filter_report.tsv
    touch ${assembly_id}.repeat_vs_tepeps.tsv ${assembly_id}.repeat_vs_prot.tsv
    touch ${assembly_id}.repeat_filter_sources.json
    """
}
