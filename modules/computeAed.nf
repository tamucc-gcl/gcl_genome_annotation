/*
========================================================================================
    COMPUTE ANNOTATION EDIT DISTANCE (AED)
========================================================================================
    Measures per-gene evidence support by aligning the same protein database
    used for BRAKER training back to the genome with miniprot, then computing
    nucleotide-level congruency between each BRAKER CDS model and the best
    overlapping protein alignment.

    AED = 1 - (SN + SP) / 2

    where:
      SN = |CDS_braker ∩ CDS_evidence| / |CDS_evidence|   (sensitivity)
      SP = |CDS_braker ∩ CDS_evidence| / |CDS_braker|     (specificity)

    A value of 0 means the gene model perfectly matches the evidence;
    1 means no evidence support at all.  Values ≤ 0.5 are considered
    well-supported; ≤ 0.3 high-quality (Holt & Yandell 2011).

    If no protein database is provided (empty file), the process emits a
    header-only AED TSV so that all downstream modules can consume the file
    without any conditional logic.

    Reference:
      Eilbeck et al. (2009) BMC Bioinformatics 10:67
      doi:10.1186/1471-2105-10-67

    Inputs:
        - assembly_id:       Sample identifier
        - masked_assembly:   Soft-masked genome assembly FASTA
        - braker_gff3:       BRAKER3 gene models GFF3
        - protein_database:  Combined protein FASTA (OrthoDB + ref proteome)
                             May be an empty file if no protein DB is available.

    Outputs:
        - aed_scores:   ${assembly_id}.aed_scores.tsv   (per-transcript AED)
        - miniprot_gff: ${assembly_id}.miniprot.gff      (protein alignments)

    Tools: miniprot, Python 3 (stdlib only)
========================================================================================
*/

process COMPUTE_AED {
    tag "${assembly_id}"
    label 'computeAed'
    publishDir "${params.outdir}/structural_annotation", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(masked_assembly), path(braker_gff3)
    path protein_database

    output:
    tuple val(assembly_id), path("${assembly_id}.aed_scores.tsv"), emit: aed_scores
    tuple val(assembly_id), path("${assembly_id}.miniprot.gff"),   emit: miniprot_gff

    script:
    """
    # ── Guard: if protein DB is empty, emit header-only stub ─────────────
    if [ ! -s "${protein_database}" ]; then
        echo "No protein database provided — emitting stub AED" >&2
        printf "gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\tprotein_identity\\tn_evidence_hits\\n" \\
            > ${assembly_id}.aed_scores.tsv
        touch ${assembly_id}.miniprot.gff
        exit 0
    fi

    # ── 1. Align proteins to genome with miniprot ────────────────────────
    miniprot \\
        -t ${task.cpus} \\
        --gff \\
        ${masked_assembly} \\
        ${protein_database} \\
        > ${assembly_id}.miniprot.gff

    # ── 2. Compute per-transcript AED ────────────────────────────────────
    python3 << 'PYEOF'
import sys
from collections import defaultdict

ASM = "${assembly_id}"
BRAKER_GFF = "${braker_gff3}"
MINIPROT_GFF = f"{ASM}.miniprot.gff"

# ── Interval utilities (0-based half-open) ────────────────────────────────

def merge_intervals(ivs):
    if not ivs:
        return []
    s = sorted(ivs)
    merged = [list(s[0])]
    for a, b in s[1:]:
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [tuple(x) for x in merged]

def iv_length(ivs):
    return sum(b - a for a, b in ivs)

def iv_intersect(a_ivs, b_ivs):
    result = []
    i = j = 0
    while i < len(a_ivs) and j < len(b_ivs):
        lo = max(a_ivs[i][0], b_ivs[j][0])
        hi = min(a_ivs[i][1], b_ivs[j][1])
        if lo < hi:
            result.append((lo, hi))
        if a_ivs[i][1] < b_ivs[j][1]:
            i += 1
        else:
            j += 1
    return result

def compute_aed(pred_ivs, evid_ivs):
    pm = merge_intervals(pred_ivs)
    em = merge_intervals(evid_ivs)
    p_len = iv_length(pm)
    e_len = iv_length(em)
    if p_len == 0 or e_len == 0:
        return 1.0
    isect = iv_intersect(pm, em)
    i_len = iv_length(isect)
    sn = i_len / e_len
    sp = i_len / p_len
    return round(1.0 - (sn + sp) / 2.0, 6)

# ── Parse GFF3 attributes ────────────────────────────────────────────────

def parse_attrs(s):
    d = {}
    for kv in s.strip().rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d

# ── Parse miniprot GFF3 ──────────────────────────────────────────────────

print("Parsing miniprot alignments...", file=sys.stderr)
mp_mrnas = {}
mp_by_scaff = defaultdict(list)

with open(MINIPROT_GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue
        scaff, src, ftype, start, end, score, strand, frame, attr_str = cols[:9]
        s0 = int(start) - 1
        e0 = int(end)
        attrs = parse_attrs(attr_str)

        if ftype == "mRNA":
            mid = attrs.get("ID", "")
            target = attrs.get("Target", "")
            prot_id = target.split()[0] if target else "."
            identity = attrs.get("Identity", ".")
            mp_mrnas[mid] = {
                "scaff": scaff, "strand": strand,
                "prot": prot_id, "ident": identity,
                "s": s0, "e": e0, "cds": [],
            }
        elif ftype == "CDS":
            parent = attrs.get("Parent", "")
            if parent in mp_mrnas:
                mp_mrnas[parent]["cds"].append((s0, e0))

for mid, info in mp_mrnas.items():
    if info["cds"]:
        info["cds"] = merge_intervals(info["cds"])
        mp_by_scaff[info["scaff"]].append((mid, info))

for scaff in mp_by_scaff:
    mp_by_scaff[scaff].sort(key=lambda x: x[1]["s"])

n_mp = sum(len(v) for v in mp_by_scaff.values())
print(f"  {n_mp} protein alignments on {len(mp_by_scaff)} scaffolds",
      file=sys.stderr)

# ── Parse BRAKER GFF3 ────────────────────────────────────────────────────

print("Parsing BRAKER gene models...", file=sys.stderr)
genes = {}
mrna_info = {}
mrna_to_gene = {}

with open(BRAKER_GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 9:
            continue
        scaff, src, ftype, start, end, score, strand, frame, attr_str = cols[:9]
        s0 = int(start) - 1
        e0 = int(end)
        attrs = parse_attrs(attr_str)

        if ftype == "gene":
            gid = attrs.get("ID", "")
            genes[gid] = {
                "scaff": scaff, "strand": strand,
                "s": s0, "e": e0, "mrnas": {},
            }
        elif ftype in ("mRNA", "transcript"):
            tid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            mrna_to_gene[tid] = parent
            info = {
                "scaff": scaff, "strand": strand,
                "s": s0, "e": e0, "cds": [],
            }
            mrna_info[tid] = info
            if parent in genes:
                genes[parent]["mrnas"][tid] = info
        elif ftype == "CDS":
            parent = attrs.get("Parent", "")
            if parent in mrna_info:
                mrna_info[parent]["cds"].append((s0, e0))

for tid in mrna_info:
    if mrna_info[tid]["cds"]:
        mrna_info[tid]["cds"] = merge_intervals(mrna_info[tid]["cds"])

n_genes = len(genes)
n_tx = sum(len(g["mrnas"]) for g in genes.values())
print(f"  {n_genes} genes, {n_tx} transcripts", file=sys.stderr)

# ── Compute AED per transcript ────────────────────────────────────────────

print("Computing AED scores...", file=sys.stderr)

results = []

for gid in sorted(genes):
    gene = genes[gid]
    scaff = gene["scaff"]
    strand = gene["strand"]
    mp_alns = mp_by_scaff.get(scaff, [])

    for tid in sorted(gene["mrnas"]):
        tinfo = gene["mrnas"][tid]
        pred_cds = tinfo["cds"]
        if not pred_cds:
            results.append((gid, tid, 1.0, ".", ".", 0))
            continue

        pred_s = pred_cds[0][0]
        pred_e = pred_cds[-1][1]

        best_aed = 1.0
        best_prot = "."
        best_ident = "."
        n_hits = 0

        for mp_id, mp_info in mp_alns:
            if mp_info["e"] <= pred_s or mp_info["s"] >= pred_e:
                continue
            if mp_info["strand"] != strand:
                continue
            if not mp_info["cds"]:
                continue

            aed = compute_aed(pred_cds, mp_info["cds"])
            n_hits += 1

            if aed < best_aed:
                best_aed = aed
                best_prot = mp_info["prot"]
                best_ident = mp_info["ident"]

        results.append((gid, tid, best_aed, best_prot, best_ident, n_hits))

# ── Write per-transcript AED table ────────────────────────────────────────

out_path = f"{ASM}.aed_scores.tsv"
with open(out_path, "w") as out:
    out.write("gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\t"
              "protein_identity\\tn_evidence_hits\\n")
    for gid, tid, aed, prot, ident, nh in results:
        out.write(f"{gid}\\t{tid}\\t{aed:.4f}\\t{prot}\\t{ident}\\t{nh}\\n")

# ── Summary statistics ────────────────────────────────────────────────────

aed_vals = [r[2] for r in results]
n_total = len(aed_vals)
if n_total == 0:
    print("WARNING: no transcripts found", file=sys.stderr)
    sys.exit(0)

n_perfect  = sum(1 for a in aed_vals if a <= 0.001)
n_lt03     = sum(1 for a in aed_vals if a < 0.3)
n_lt05     = sum(1 for a in aed_vals if a < 0.5)
n_no_evid  = sum(1 for a in aed_vals if a >= 0.999)
mean_aed   = sum(aed_vals) / n_total

gene_aeds = {}
for gid, tid, aed, *_ in results:
    if gid not in gene_aeds or aed < gene_aeds[gid]:
        gene_aeds[gid] = aed

gv = list(gene_aeds.values())
ng = len(gv)
g03 = sum(1 for a in gv if a < 0.3)
g05 = sum(1 for a in gv if a < 0.5)

pct = lambda n, d: f"{100*n/d:.1f}%" if d else "N/A"

print(f"\\nAED Summary ({ASM}):", file=sys.stderr)
print(f"  Transcripts:     {n_total}", file=sys.stderr)
print(f"  AED = 0:         {n_perfect} ({pct(n_perfect, n_total)})", file=sys.stderr)
print(f"  AED < 0.3:       {n_lt03} ({pct(n_lt03, n_total)})", file=sys.stderr)
print(f"  AED < 0.5:       {n_lt05} ({pct(n_lt05, n_total)})", file=sys.stderr)
print(f"  AED = 1 (no ev): {n_no_evid} ({pct(n_no_evid, n_total)})", file=sys.stderr)
print(f"  Mean AED:        {mean_aed:.4f}", file=sys.stderr)
print(f"  Genes:", file=sys.stderr)
print(f"    AED < 0.3:     {g03}/{ng} ({pct(g03, ng)})", file=sys.stderr)
print(f"    AED < 0.5:     {g05}/{ng} ({pct(g05, ng)})", file=sys.stderr)
PYEOF
    """

    stub:
    """
    touch ${assembly_id}.miniprot.gff
    printf "gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\tprotein_identity\\tn_evidence_hits\\n" \\
        > ${assembly_id}.aed_scores.tsv
    """
}


process AED_STUB {
    tag "${assembly_id}"
    label 'mergeFiles'

    input:
    tuple val(assembly_id), path(braker_gff3)

    output:
    tuple val(assembly_id), path("${assembly_id}.aed_scores.tsv"), emit: aed_scores

    script:
    """
    printf "gene_id\\ttranscript_id\\taed\\tbest_protein_hit\\tprotein_identity\\tn_evidence_hits\\n" \\
        > ${assembly_id}.aed_scores.tsv
    """
}