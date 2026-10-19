#!/usr/bin/env python3
"""Markdown/status reporting only; changing this helper does not rerun native searches."""
import argparse
import base64
import json
from pathlib import Path
from genome_tracks import TRACKS, dictionary
from reference_repeats import dump


def report(aid, settings, identity, statuses, outdir):
    data, _ = dictionary(identity)
    root = Path(outdir); root.mkdir(parents=True, exist_ok=True)
    states = [json.loads(Path(p).read_text()) for p in statuses]
    expected = {t for t in TRACKS if settings['run_' + t]}
    if len(states) != len(expected) or {s['track'] for s in states} != expected:
        raise ValueError('Missing/duplicate/unexpected genome track results')
    if any(s['status'] != 'ready' or s['assembly_id'] != aid or s['reference_identity'] != data['reference_identity'] for s in states):
        raise ValueError('Genome track/reference identity mismatch')
    by_track = {s['track']: s for s in states}
    lines = [f'# Genome annotation tracks: {aid}\n\n', 'Status: **ready**\n\n',
             '| Track | Status | Features | Validated GFF3 |\n|---|---|---:|---|\n']
    for track in TRACKS:
        state = by_track.get(track)
        lines.append(f'| {track} | {"ready" if state else "disabled"} | {state["features"] if state else "—"} | ' +
                     (f'[{track}]({aid}.{track}.gff3)' if state else '—') + ' |\n')
    lines += ['\nZero features is a valid completed search. Feature counts include child features where present.\n\n',
              'tRNAscan and Rfam are separate evidence tracks and can overlap; counts are not a deduplicated ncRNA gene total.\n\n',
              'Telomeric regions are window-based motif predictions. Centromeric satellite and higher-order repeat candidates do not establish functional centromeres.\n\n',
              'Raw tool results: [ncRNA](../../ncRNA/), [telomeres](../../telomeres/), [centromere](../../centromere/). Coding-gene identifiers and structures are unchanged.\n']
    (root/'README.md').write_text(''.join(lines), encoding='utf-8')
    state = dict(schema_version=1, stage='genome_tracks', assembly_id=aid, status='ready',
                 reference_identity=data['reference_identity'], settings=settings,
                 tracks={t: by_track.get(t, dict(status='disabled')) for t in TRACKS})
    dump(root/'genome_tracks_status.json', state)
    return state



if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('aid','settings','identity','outdir'):p.add_argument('--'+name,required=True)
    p.add_argument('--statuses',nargs='*',default=[])
    a=vars(p.parse_args());a['settings']=json.loads(base64.b64decode(a['settings']))
    report(**a)
