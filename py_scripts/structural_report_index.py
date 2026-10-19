"""Write the Markdown assembly index without embedded shell/Python quoting."""
import argparse
import base64
import json
from pathlib import Path


def write_index(reports,expected,out):
    rows=[json.loads((Path(p)/'structural_qc_status.json').read_text()) for p in reports]
    if sorted(r['assembly_id'] for r in rows)!=sorted(expected):
        raise ValueError('Missing/duplicate assembly structural reports')
    if any(r.get('status')!='ready' for r in rows):raise ValueError('Structural report is not ready')
    lines=['# Structural annotation QC','','| Assembly | Before genes | After genes | Before transcripts | After transcripts |',
           '|---|---:|---:|---:|---:|']
    for row in sorted(rows,key=lambda r:r['assembly_id']):
        aid=row['assembly_id'];before,after=row['model_sets']
        lines.append(f"| [{aid}]({aid}/README.md) | {before['genes']} | {after['genes']} | {before['transcripts']} | {after['transcripts']} |")
    Path(out).write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--expected',required=True)
    p.add_argument('--reports',nargs='+',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    write_index(a.reports,json.loads(base64.b64decode(a.expected,validate=True)),a.out)
