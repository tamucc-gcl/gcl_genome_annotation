#!/usr/bin/env python3
"""Parsers retained from MERGE_FUNCTIONAL_ANNOTATIONS; joins accepted assessment."""
import argparse
import json
from pathlib import Path
from functional_inputs import TOOLS, read_table
from prepare_protein_evidence import checksum

FUNCTIONAL_FIELDS = '''ipr_signatures ipr_accessions ipr_descriptions ipr_go
eggnog_OGs eggnog_name eggnog_desc eggnog_GO eggnog_KEGG_ko eggnog_KEGG_path eggnog_CAZy eggnog_PFAMs
diamond_hit diamond_desc diamond_pident diamond_evalue
kofam_KO kofam_def kofam_score kofam_evalue kofam_sig
deepsig_type deepsig_start deepsig_end deepsig_score deepsig_attr'''.split()


def parse_results(args):
    import csv, re, sys, statistics
    from collections import defaultdict, OrderedDict

    # ── Storage: final renamed protein_id → dict of annotation fields ─────────
    annot = defaultdict(lambda: {})

    def parse_gff_attrs(attr_str):
        attrs = {}
        if not attr_str or attr_str == ".":
            return attrs
        for tok in attr_str.split(";"):
            if "=" in tok:
                k, v = tok.split("=", 1)
                attrs[k.strip()] = v.strip()
        return attrs

    # ── 1. InterProScan (headerless TSV, 11-15 columns) ─────────────────────
    #   Col layout (0-indexed):
    #     0  protein_id
    #     3  analysis (e.g. Pfam, PANTHER)
    #     4  signature_accession
    #     5  signature_description
    #     11 IPR accession  (may be absent)
    #     12 IPR description (may be absent)
    #     13 GO terms        (may be absent)
    with open(args.interproscan_tsv) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11: raise ValueError("Malformed InterProScan TSV")
            pid = cols[0]
            rec = annot[pid]

            sig = cols[4] if len(cols) > 4 else "-"
            sig_desc = cols[5] if len(cols) > 5 else "-"
            ipr = cols[11] if len(cols) > 11 and cols[11] else "-"
            ipr_desc = cols[12] if len(cols) > 12 and cols[12] else "-"
            go = cols[13] if len(cols) > 13 and cols[13] else "-"

            entry = f"{sig}:{sig_desc}"
            rec.setdefault("ipr_signatures", []).append(entry)
            if ipr != "-":
                rec.setdefault("ipr_accessions", set()).add(ipr)
            if ipr_desc != "-":
                rec.setdefault("ipr_descriptions", set()).add(ipr_desc)
            if go != "-":
                for g in go.split("|"):
                    rec.setdefault("ipr_go", set()).add(g.strip())

    # ── 2. eggNOG (header starts with #, data lines are tab-separated) ──────
    with open(args.eggnog_tsv) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 21: raise ValueError("Malformed eggNOG TSV")
            pid = cols[0]
            rec = annot[pid]
            rec["eggnog_OGs"]       = cols[4]  if len(cols) > 4  and cols[4]  else "-"
            rec["eggnog_desc"]      = cols[7]  if len(cols) > 7  and cols[7]  else "-"
            rec["eggnog_name"]      = cols[8]  if len(cols) > 8  and cols[8]  else "-"
            rec["eggnog_GO"]        = cols[9]  if len(cols) > 9  and cols[9]  else "-"
            rec["eggnog_KEGG_ko"]   = cols[11] if len(cols) > 11 and cols[11] else "-"
            rec["eggnog_KEGG_path"] = cols[12] if len(cols) > 12 and cols[12] else "-"
            rec["eggnog_CAZy"]      = cols[15] if len(cols) > 15 and cols[15] else "-"
            rec["eggnog_PFAMs"]     = cols[20] if len(cols) > 20 and cols[20] else "-"

    # ── 3. DIAMOND (first line is header) ────────────────────────────────────
    with open(args.diamond_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row.get("qseqid", "").strip()
            if not pid:
                continue
            rec = annot[pid]
            if "diamond_hit" not in rec:
                rec["diamond_hit"]    = row.get("sseqid", "-")
                rec["diamond_desc"]   = row.get("stitle", "-")
                rec["diamond_pident"] = row.get("pident", "-")
                rec["diamond_evalue"] = row.get("evalue", "-")

    # ── 4. KofamScan (first line is header, detail-tsv format) ───────────────
    with open(args.kofamscan_tsv) as fh:
        for i, line in enumerate(fh):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                raise ValueError("Malformed KofamScan detail-tsv")
            sig_marker = cols[0].strip()
            pid = cols[1].strip()
            ko = cols[2].strip() if cols[2].strip() else "-"
            score = cols[4].strip() if len(cols) > 4 else "-"
            evalue = cols[5].strip() if len(cols) > 5 else "-"
            ko_def = cols[6].strip() if len(cols) > 6 else "-"

            if sig_marker != "*":
                continue
            rec = annot[pid]
            if "kofam_KO" not in rec or (float(evalue), -float(score), ko) < (float(rec["kofam_evalue"]), -float(rec["kofam_score"]), rec["kofam_KO"]):
                rec["kofam_KO"]     = ko
                rec["kofam_score"]  = score
                rec["kofam_evalue"] = evalue
                rec["kofam_def"]    = ko_def
                rec["kofam_sig"]    = sig_marker if sig_marker else "-"

    # ── 5. DeepSig (GFF3 format) ────────────────────────────────────────────
    with open(args.deepsig_gff3) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise ValueError("Malformed DeepSig GFF")
            pid = cols[0]
            feat_type = cols[2] if len(cols) > 2 else "-"
            start = cols[3] if len(cols) > 3 else "-"
            end = cols[4] if len(cols) > 4 else "-"
            score = cols[5] if len(cols) > 5 else "-"
            attributes = cols[8] if len(cols) > 8 else "-"

            if feat_type.lower().replace(" ", "_") not in ("signal_peptide", "signalpeptide", "transmembrane"):
                continue
            rec = annot[pid]
            if "deepsig_type" not in rec:
                rec["deepsig_type"]  = feat_type
                rec["deepsig_start"] = start
                rec["deepsig_end"]   = end
                rec["deepsig_score"] = score
                rec["deepsig_attr"]  = attributes
    return annot


def merge(args):
    import csv
    from collections import Counter
    fields, assessment = read_table(args.assessment, 'transcript_id')
    _, eligible = read_table(args.eligibility, 'protein_id')
    if set(assessment) != set(eligible):
        raise ValueError('Eligibility/assessment IDs differ')
    import base64
    settings = json.loads(base64.b64decode(args.settings))
    annot = parse_results(args)
    extra = set(annot)-set(assessment)
    if extra:
        raise ValueError('Unknown functional protein IDs: '+', '.join(sorted(extra)[:10]))
    if any(eligible[t]['eligible'] != 'true' for t in annot):
        raise ValueError('Functional result for an ineligible protein')
    prefixes = dict(interproscan='ipr_', eggnog='eggnog_', diamond='diamond_',
                    kofamscan='kofam_', deepsig='deepsig_')
    counts = {tool: Counter() for tool in TOOLS}
    rows = []
    for tid, metric in sorted(assessment.items()):
        row = dict(protein_id=tid, **metric)
        rec = annot.get(tid, {})
        for key, value in rec.items():
            row[key] = '; '.join(sorted(set(value))) if isinstance(value, (set,list)) else value
        for tool in TOOLS:
            hit = any(k.startswith(prefixes[tool]) and v not in ('', '.', '-') for k,v in rec.items())
            enabled = settings.get('run_'+tool, False)
            if hit and not enabled:
                raise ValueError('Results supplied for disabled '+tool)
            state = ('disabled' if not enabled else 'ineligible' if eligible[tid]['eligible']!='true'
                     else 'hit' if hit else 'no_hit')
            row[tool+'_status'] = state
            counts[tool][state] += 1
        row.update(search_eligible=eligible[tid]['eligible'], search_ineligible_reason=eligible[tid]['reason'])
        rows.append(row)
    # OMArk classified only the selected representative proteins, not every isoform.
    classes = {}; category = None
    if args.omark:
        for line in Path(args.omark).read_text().splitlines():
            line = line.strip()
            if line.startswith('>'): category = line[1:]
            elif line and not line.startswith('#'):
                if line not in assessment or line in classes or category is None:
                    raise ValueError('Invalid/duplicate selected OMArk ID: '+line)
                classes[line] = category
    for row in rows:
        row['omark_class'] = classes.get(row['protein_id'], '.')
        row['omark_status'] = ('classified' if row['protein_id'] in classes else
                               'not_classified' if args.omark else 'disabled')
    columns = ['protein_id']+fields+FUNCTIONAL_FIELDS+sorted(
        set().union(*(r.keys() for r in rows))-set(fields)-set(FUNCTIONAL_FIELDS)-{'protein_id'})
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    target = out/f'{args.assembly_id}.functional_annotations.tsv'
    with target.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t'); writer.writeheader()
        writer.writerows({k:r.get(k,'.') for k in columns} for r in rows)
    result = dict(status='ready', stage='functional_merge', assembly_id=args.assembly_id,
                  transcripts=len(rows), genes=len({r['gene_id'] for r in rows}), tools=counts,
                  source_assessment_sha256=checksum(args.assessment), table_sha256=checksum(target),
                  settings=settings, sources={k:checksum(getattr(args,k)) for k in
                      ('interproscan_tsv','eggnog_tsv','diamond_tsv','kofamscan_tsv','deepsig_gff3')})
    (out/'functional_merge_status.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('assembly-id','assessment','eligibility','settings','outdir',
                'interproscan-tsv','eggnog-tsv','diamond-tsv','kofamscan-tsv','deepsig-gff3'):
        p.add_argument('--'+key, required=True)
    p.add_argument('--omark')
    merge(p.parse_args())
