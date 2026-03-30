#!/usr/bin/env python3
"""
build_coord_map.py
==================
Build coord_map.tsv for REPEAT_MASK_BATCH by parsing actual FASTA headers
from the batch fasta produced by seqkit subseq --bed.

seqkit subseq --bed outputs sequences in genome fasta order, NOT BED file
order, so the coord_map must be derived from the actual header content rather
than assuming BED line N == fasta sequence N.

Header format produced by seqkit subseq --bed:
    >scaffold_name_start-end:.
where start and end are 1-based coordinates.

Output coord_map.tsv columns:
    short_id  scaffold_name  1based_start  end  scaffold_length

Usage:
    python3 build_coord_map.py <batch_fasta> <seq_lengths_tsv> <output_tsv>
"""
import re
import sys

batch_fasta   = sys.argv[1]
seq_lengths   = sys.argv[2]
output_tsv    = sys.argv[3]

# Load scaffold lengths
slen = {}
with open(seq_lengths) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            slen[parts[0]] = int(parts[1])

# Parse fasta headers and write coord_map
n = 0
with open(batch_fasta) as fh, open(output_tsv, "w") as out:
    for line in fh:
        if not line.startswith(">"):
            continue
        short = f"seq_{n:04d}"
        n += 1
        # Strip > and anything after whitespace, then strip :. suffix
        hdr = line[1:].split()[0]
        hdr = hdr.split(":")[0]
        # Match scaffold_name_start-end pattern
        m = re.match(r'^(.+)_([0-9]+)-([0-9]+)$', hdr)
        if m:
            scaff     = m.group(1)
            start1    = int(m.group(2))
            end_coord = int(m.group(3))
        else:
            scaff     = hdr
            start1    = 1
            end_coord = slen.get(hdr, 0)
        scaffold_len = slen.get(scaff, 0)
        out.write(f"{short}\t{scaff}\t{start1}\t{end_coord}\t{scaffold_len}\n")

print(f"coord_map.tsv: {n} sequences mapped", file=sys.stderr)