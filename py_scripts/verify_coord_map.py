#!/usr/bin/env python3
"""
verify_coord_map.py
===================
Verify that each sequence in the adjusted fasta has the correct length
by comparing actual sequence length against the coordinates encoded in
the header (scaffold__start__end).

Prints a summary to stderr. Any MISMATCH line indicates a coord_map
ordering error that needs investigation.

Usage:
    python3 verify_coord_map.py <adjusted_fasta>
"""
import sys

adjusted_fa = sys.argv[1]

print(f"=== COORD_MAP VERIFICATION: {adjusted_fa} ===", file=sys.stderr)

mismatches = 0
checked = 0
name = None
seq_len = 0

def check(name, seq_len):
    parts = name.rsplit("__", 2)
    if len(parts) == 3:
        try:
            expected = int(parts[2]) - int(parts[1]) + 1
            status = "OK" if seq_len == expected else f"MISMATCH(got={seq_len},expected={expected})"
            return status, seq_len == expected
        except ValueError:
            pass
    return "MISMATCH(missing_range)", False

with open(adjusted_fa) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                status, ok = check(name, seq_len)
                print(f"  {name}\tlength={seq_len}\t{status}", file=sys.stderr)
                checked += 1
                if not ok:
                    mismatches += 1
            name = line[1:]
            seq_len = 0
        else:
            seq_len += len(line)

if name is not None:
    status, ok = check(name, seq_len)
    print(f"  {name}\tlength={seq_len}\t{status}", file=sys.stderr)
    checked += 1
    if not ok:
        mismatches += 1

print(f"=== {checked} sequences checked, {mismatches} mismatches ===", file=sys.stderr)
if mismatches > 0:
    sys.exit(1)
