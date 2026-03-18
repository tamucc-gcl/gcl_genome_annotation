/*
========================================================================================
    EXTRACT SEQUENCES MODULE
========================================================================================
    Extract protein, CDS, transcript, and tRNA FASTA sequences from the
    final annotated GFF3 and genome assembly.

    This produces the canonical sequence files with pipeline-namespaced IDs
    ({assembly_id}_g{N}.t{M} for protein-coding, {assembly_id}_tRNA{N} for
    tRNAs) and functional product names in FASTA headers.

    Protein-coding sequences are extracted with gffread from BRAKER3-sourced
    features. tRNA sequences are extracted via a Python script that handles
    intron splicing and builds descriptive headers with isotype, anticodon,
    and genomic coordinates.

    Inputs:
        - assembly_id:      Sample identifier
        - annotated_gff3:   Merged annotated GFF3 (from ANNOTATE_GFF3)
        - genome_fasta:     Genome assembly FASTA (filtered, unmasked)

    Outputs:
        - proteins:     ${assembly_id}.proteins.fa     (translated protein sequences)
        - cds:          ${assembly_id}.cds.fa           (CDS nucleotide sequences)
        - transcripts:  ${assembly_id}.transcripts.fa   (spliced mRNA sequences)
        - trna:         ${assembly_id}.tRNA.fa           (spliced tRNA sequences)

    Tools: gffread (conda: bioconda::gffread), Python 3 (stdlib only)
========================================================================================
*/

process EXTRACT_SEQUENCES {
    tag "${assembly_id}"
    label 'extract_sequences'
    publishDir "${params.outdir}/sequences", mode: params.publish_dir_mode

    input:
    tuple val(assembly_id), path(annotated_gff3), path(genome_fasta)

    output:
    tuple val(assembly_id), path("${assembly_id}.proteins.fa"),    emit: proteins
    tuple val(assembly_id), path("${assembly_id}.cds.fa"),         emit: cds
    tuple val(assembly_id), path("${assembly_id}.transcripts.fa"), emit: transcripts
    tuple val(assembly_id), path("${assembly_id}.tRNA.fa"),        emit: trna

    script:
    """
    set -euo pipefail

    # ── 1. Extract protein-coding features from annotated GFF3 ────────────────
    #   Filter to BRAKER3-sourced lines only (gene, mRNA, exon, CDS, etc.)
    awk '/^#/ || \$2 == "BRAKER3"' ${annotated_gff3} > coding_raw.gff3

    # ── 1b. Validate coordinates against scaffold lengths ─────────────────────
    #   BRAKER can occasionally produce CDS features extending past scaffold
    #   ends. gffread aborts on these, so we identify and remove the offending
    #   gene models (entire gene tree) before extraction.
    seqkit fx2tab -n -l ${genome_fasta} > scaffold_lengths.tsv

    python3 - <<'VALIDATE'
import sys

# Load scaffold lengths
slen = {}
with open("scaffold_lengths.tsv") as fh:
    for line in fh:
        parts = line.strip().split("\t")
        slen[parts[0]] = int(parts[1])

# Pass 1: find all features with out-of-bounds coordinates, collect their
# gene-level parent IDs so we can remove entire gene trees
bad_genes = set()
features = []

with open("coding_raw.gff3") as fh:
    for line in fh:
        if line.startswith("#"):
            features.append(("comment", line, None, None))
            continue
        line_s = line.strip()
        if not line_s:
            continue
        cols = line_s.split("\t")
        if len(cols) != 9:
            features.append(("comment", line, None, None))
            continue

        seqid = cols[0]
        end   = int(cols[4])
        attrs = {}
        for token in cols[8].split(";"):
            if "=" in token:
                k, v = token.split("=", 1)
                attrs[k.strip()] = v.strip()

        scaffold_len = slen.get(seqid, 0)
        oob = end > scaffold_len

        features.append(("feature", line, attrs, oob))

        if oob:
            ftype = cols[2]
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            # Walk up to gene level
            if ftype == "gene":
                bad_genes.add(fid)
            elif ftype in ("mRNA", "transcript"):
                bad_genes.add(parent)
                bad_genes.add(fid)  # also mark mRNA so children get caught
            else:
                # CDS/exon — parent is mRNA, grandparent is gene
                bad_genes.add(parent)

# Pass 2: find mRNA parents of bad mRNAs → bad genes
# Re-scan to resolve mRNA→gene mappings
mrna_to_gene = {}
for tag, line, attrs, oob in features:
    if tag != "feature" or attrs is None:
        continue
    cols = line.strip().split("\t")
    if cols[2] in ("mRNA", "transcript"):
        mrna_to_gene[attrs.get("ID", "")] = attrs.get("Parent", "")

# Expand bad_genes: if a bad mRNA ID is in the set, add its gene parent
expanded = set()
for bg in bad_genes:
    expanded.add(bg)
    if bg in mrna_to_gene:
        expanded.add(mrna_to_gene[bg])
bad_genes = expanded

n_removed = 0
with open("coding.gff3", "w") as out:
    for tag, line, attrs, oob in features:
        if tag == "comment":
            out.write(line)
            continue
        fid = attrs.get("ID", "") if attrs else ""
        parent = attrs.get("Parent", "") if attrs else ""
        # Skip if this feature or its parent belongs to a bad gene tree
        if fid in bad_genes or parent in bad_genes:
            n_removed += 1
            continue
        out.write(line)

if n_removed > 0:
    print(f"WARNING: Removed {n_removed} GFF3 features from "
          f"{len(bad_genes)} gene models with out-of-bounds coordinates",
          file=sys.stderr)
VALIDATE

    # ── 2. Protein sequences (translated from CDS) ───────────────────────────
    gffread coding.gff3 -g ${genome_fasta} -y ${assembly_id}.proteins.fa -F -S

    # ── 3. CDS nucleotide sequences ──────────────────────────────────────────
    gffread coding.gff3 -g ${genome_fasta} -x ${assembly_id}.cds.fa -F

    # ── 4. Spliced transcript (mRNA) sequences ───────────────────────────────
    gffread coding.gff3 -g ${genome_fasta} -w ${assembly_id}.transcripts.fa -F

    # ── 5. tRNA sequences (Python — handles intron splicing + rich headers) ──
    python3 - <<'PYEOF'
import sys
from collections import OrderedDict

ASM = "${assembly_id}"

# ── Parse genome FASTA into memory ───────────────────────────────────────
genome = {}
current_id = None
current_seq = []

with open("${genome_fasta}") as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if current_id is not None:
                genome[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.upper())
    if current_id is not None:
        genome[current_id] = "".join(current_seq)

# ── Parse annotated GFF3 for tRNAscan-SE features ───────────────────────
#   Build hierarchy: gene → tRNA transcript → exon(s)

def parse_attrs(attr_str):
    d = OrderedDict()
    if not attr_str or attr_str == ".":
        return d
    for token in attr_str.split(";"):
        token = token.strip()
        if "=" in token:
            k, v = token.split("=", 1)
            d[k.strip()] = v.strip()
    return d

trna_genes = {}        # gene_id → {seqid, start, end, strand, attrs}
trna_transcripts = {}  # trna_id → {parent, seqid, start, end, strand, attrs}
trna_exons = {}        # trna_id → [(start, end), ...]

with open("${annotated_gff3}") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.strip().split("\\t")
        if len(cols) != 9 or cols[1] != "tRNAscan-SE":
            continue

        seqid  = cols[0]
        ftype  = cols[2]
        start  = int(cols[3])   # 1-based inclusive
        end    = int(cols[4])   # 1-based inclusive
        strand = cols[6]
        attrs  = parse_attrs(cols[8])

        if ftype == "gene":
            trna_genes[attrs.get("ID", "")] = {
                "seqid": seqid, "start": start, "end": end,
                "strand": strand, "attrs": attrs,
            }
        elif ftype == "tRNA":
            tid = attrs.get("ID", "")
            trna_transcripts[tid] = {
                "parent": attrs.get("Parent", ""),
                "seqid": seqid, "start": start, "end": end,
                "strand": strand, "attrs": attrs,
            }
            trna_exons.setdefault(tid, [])
        elif ftype == "exon":
            parent = attrs.get("Parent", "")
            if parent in trna_exons:
                trna_exons[parent].append((start, end))

# ── Reverse complement helper ────────────────────────────────────────────
COMP = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")

def revcomp(seq):
    return seq.translate(COMP)[::-1]

# ── Extract and write tRNA sequences ─────────────────────────────────────
n_written = 0
n_intron = 0

with open(f"{ASM}.tRNA.fa", "w") as out:
    # Sort by genomic position for deterministic output
    sorted_tids = sorted(
        trna_transcripts.keys(),
        key=lambda t: (
            trna_transcripts[t]["seqid"],
            trna_transcripts[t]["start"],
        ),
    )

    for tid in sorted_tids:
        info = trna_transcripts[tid]
        parent_id = info["parent"]
        gene_info = trna_genes.get(parent_id, {})

        seqid  = info["seqid"]
        strand = info["strand"]

        # Exon coordinates (sorted by genomic position)
        exons = sorted(trna_exons.get(tid, []), key=lambda x: x[0])

        # Fall back to tRNA feature span if no exon children
        if not exons:
            exons = [(info["start"], info["end"])]

        has_intron = len(exons) > 1
        if has_intron:
            n_intron += 1

        # Extract and concatenate exon sequences (GFF3 = 1-based inclusive)
        scaffold_seq = genome.get(seqid, "")
        if not scaffold_seq:
            print(f"WARNING: scaffold {seqid} not found in genome, "
                  f"skipping {tid}", file=sys.stderr)
            continue

        seq_parts = []
        for ex_start, ex_end in exons:
            seq_parts.append(scaffold_seq[ex_start - 1 : ex_end])

        seq = "".join(seq_parts)

        if strand == "-":
            seq = revcomp(seq)

        if not seq:
            continue

        # ── Build descriptive FASTA header ────────────────────────────
        #   Format: >ID Name(anticodon) scaffold:start-end(strand) [exons=N]
        name      = info["attrs"].get("Name", gene_info.get("attrs", {}).get("Name", "tRNA"))
        anticodon = info["attrs"].get("anticodon", "")
        isotype   = info["attrs"].get("isotype", "")

        gene_start = gene_info.get("start", info["start"])
        gene_end   = gene_info.get("end", info["end"])

        # Prefer isotype for display if Name is generic
        display_name = name
        if anticodon:
            display_name += f"({anticodon})"

        desc = f"{display_name} {seqid}:{gene_start}-{gene_end}({strand})"

        if has_intron:
            desc += f" exons={len(exons)}"

        # Include isotype score if present
        iso_score = info["attrs"].get("isotype_score", "")
        if iso_score:
            desc += f" isotype_score={iso_score}"

        out.write(f">{tid} {desc}\\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i : i + 60] + "\\n")

        n_written += 1

# ── Summary stats ────────────────────────────────────────────────────────
print(
    f"Extracted {n_written} tRNA sequences "
    f"({n_intron} with introns)",
    file=sys.stderr,
)
PYEOF

    # ── 6. Summary ───────────────────────────────────────────────────────────
    n_prot=\$(grep -c '^>' ${assembly_id}.proteins.fa)
    n_cds=\$(grep -c '^>' ${assembly_id}.cds.fa)
    n_tx=\$(grep -c '^>' ${assembly_id}.transcripts.fa)
    n_trna=\$(grep -c '^>' ${assembly_id}.tRNA.fa)
    echo "Extracted sequences for ${assembly_id}: \${n_prot} proteins, \${n_cds} CDS, \${n_tx} transcripts, \${n_trna} tRNAs" >&2
    """

    stub:
    """
    touch ${assembly_id}.proteins.fa
    touch ${assembly_id}.cds.fa
    touch ${assembly_id}.transcripts.fa
    touch ${assembly_id}.tRNA.fa
    """
}