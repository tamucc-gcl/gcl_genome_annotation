#!/usr/bin/env python3
"""Package accepted annotation outputs; perform no searches or model selection."""
import argparse
import base64
from collections import Counter
import csv
import gc
import hashlib
import json
from pathlib import Path
import re
import shutil
from urllib.parse import quote, unquote

from reference_repeats import fasta, sha, dump
from model_metrics import sequences, revcomp
from evidence_aed import models
from functional_report import geometry


TRACKS = ('trnascan', 'infernal', 'telomere', 'centromere')


def accepted(path, aid, ref=None, stage=None):
    state = json.loads(Path(path).read_text())
    if state.get('status') != 'ready' or state.get('assembly_id') != aid:
        raise ValueError('Unaccepted annotation input: '+str(path))
    if ref is not None and state.get('reference_identity') != ref:
        raise ValueError('Reference identity mismatch: '+str(path))
    if stage and state.get('stage') != stage:
        raise ValueError('Unexpected annotation stage: '+str(path))
    return state


def gff_records(path, lengths):
    rows = []; nodes = {}; bounds = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'): continue
            c = line.rstrip('\r\n').split('\t')
            if len(c) != 9: raise ValueError('Expected GFF3, not GTF: '+str(path))
            a = dict(v.split('=', 1) for v in c[8].strip(';').split(';') if '=' in v)
            if c[0] not in lengths or not 1 <= int(c[3]) <= int(c[4]) <= lengths[c[0]] or c[6] not in ('+', '-', '.', '?'):
                raise ValueError('GFF3 coordinates/strand differ from reference: '+str(path))
            fid = a.get('ID')
            if not fid: raise ValueError('GFF3 feature has no ID: '+str(path))
            signature = (c[0], c[2], c[6], a.get('Parent', ''), a.get('Derives_from', ''))
            if fid in nodes and nodes[fid] != signature: raise ValueError('Inconsistent multipart GFF3 ID: '+fid)
            nodes[fid] = signature
            low, high = bounds.get(fid, (int(c[3]), int(c[4])))
            bounds[fid] = (min(low, int(c[3])), max(high, int(c[4])))
            rows.append((c, a))
    edges = {fid: tuple(p for value in (s[3], s[4]) for p in value.split(',') if p) for fid, s in nodes.items()}
    if any(p not in nodes for parents in edges.values() for p in parents): raise ValueError('Dangling GFF3 feature relationship')
    for fid, signature in nodes.items():
        for parent in signature[3].split(','):
            if not parent: continue
            p = nodes[parent]
            if p[0] != signature[0] or (p[2] in ('+', '-') and signature[2] in ('+', '-') and p[2] != signature[2]) or not bounds[parent][0] <= bounds[fid][0] <= bounds[fid][1] <= bounds[parent][1]:
                raise ValueError('GFF3 child differs from parent reference, strand or bounds: '+fid)
    visited = set(); visiting = set()
    def visit(fid):
        if fid in visited: return
        if fid in visiting: raise ValueError('Cyclic GFF3 feature relationship')
        visiting.add(fid)
        for parent in edges[fid]: visit(parent)
        visiting.remove(fid); visited.add(fid)
    for fid in nodes: visit(fid)
    return rows


def write_gff(path, rows, lengths):
    # Parents precede their children, including shared and multipart features.
    groups = {}; order = []
    for c, a in rows:
        if a['ID'] not in groups: groups[a['ID']] = []; order.append(a['ID'])
        groups[a['ID']].append((c, a))
    emitted = set()
    with Path(path).open('w') as handle:
        handle.write('##gff-version 3\n')
        for seq, length in lengths.items(): handle.write(f'##sequence-region {seq} 1 {length}\n')
        def emit(fid):
            if fid in emitted: return
            for key in ('Parent', 'Derives_from'):
                for parent in groups[fid][0][1].get(key, '').split(','):
                    if parent: emit(parent)
            for c, a in groups[fid]:
                handle.write('\t'.join(c[:8]+[';'.join(k+'='+v for k, v in a.items())])+'\n')
            emitted.add(fid)
        for fid in order: emit(fid)


def trna_models(rows):
    # Reuse the original extractor's exon-splicing/span-fallback behavior, but
    # use explicit tRNA IDs/parents rather than a particular source or gene tree.
    targets = {}
    for c, a in rows:
        if c[2] != 'tRNA': continue
        if a['ID'] in targets: raise ValueError('Duplicate tRNA sequence ID: '+a['ID'])
        targets[a['ID']] = dict(seq=c[0], strand=c[6], start=int(c[3])-1, end=int(c[4]), exons=[])
    for c, a in rows:
        if c[2] != 'exon': continue
        for parent in a.get('Parent', '').split(','):
            if parent in targets: targets[parent]['exons'].append((int(c[3])-1, int(c[4])))
    for m in targets.values():
        if not m['exons']: m['exons'] = [(m['start'], m['end'])]
    return targets


def extract_spliced(genome, groups):
    """Extract final-ID exon sequences in one reference pass (one contig in memory)."""
    wanted = {}; results = {}
    for label, targets, destination in groups:
        results[label] = {}
        for tid, model in targets.items():
            blocks = sorted(set(model['exons']))
            if model['strand'] not in ('+', '-') or not blocks:
                raise ValueError('Sequence export needs resolved strand/exons: '+tid)
            if any(s < model['start'] or e > model['end'] or s >= e for s,e in blocks) or any(l[1] > r[0] for l,r in zip(blocks,blocks[1:])):
                raise ValueError('Invalid/overlapping sequence export exons: '+tid)
            wanted.setdefault(model['seq'], []).append((label, tid, model['strand'], blocks))
    for header, sequence in fasta(genome):
        for label, tid, strand, blocks in wanted.pop(header.split()[0], []):
            if blocks[-1][1] > len(sequence): raise ValueError('Sequence export outside reference: '+tid)
            spliced = ''.join(sequence[s:e] for s,e in blocks).upper()
            results[label][tid] = revcomp(spliced) if strand == '-' else spliced
    if wanted: raise ValueError('Sequence export contig missing from reference')
    counts = {}
    for label, targets, destination in groups:
        if set(results[label]) != set(targets): raise ValueError('Sequence export IDs differ: '+label)
        with Path(destination).open('w') as handle:
            for tid, sequence in sorted(results[label].items()):
                handle.write('>'+tid+'\n')
                for start in range(0,len(sequence),60): handle.write(sequence[start:start+60]+'\n')
        counts[label] = len(targets)
    return counts


def package(request, outdir):
    aid = request['aid']
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', aid): raise ValueError('Invalid release assembly ID')
    out = Path(outdir)
    if out.exists() and any(out.iterdir()): raise ValueError('Release output must be empty; refusing stale files')
    out.mkdir(parents=True, exist_ok=True)
    identified = Path(request['identified']); identity = json.loads(Path(request['identity']).read_text())
    ref = identity['reference_identity']; lengths = {c['id']: c['length'] for c in identity['contigs']}
    annotation = accepted(identified/'annotation_status.json', aid, ref, 'identified_annotation')
    for name, digest in annotation['outputs'].items():
        if Path(name).name != name or sha(identified/name) != digest: raise ValueError('Identified output changed: '+name)
    for genome in (request['genome'], request['masked']):
        observed = [dict(id=h.split()[0], length=len(seq), uppercase_sha256=hashlib.sha256(seq.upper().encode('ascii')).hexdigest()) for h, seq in fasta(genome)]
        if observed != identity['contigs']: raise ValueError('Release FASTA differs from reference identity')
    repeats = accepted(request['repeat_status'], aid, ref, 'reference_repeats')
    if repeats['masked_sha256'] != sha(request['masked']): raise ValueError('Accepted masked FASTA changed')
    track_state = accepted(request['track_status'], aid, ref, 'genome_tracks')
    track_paths = {Path(p).name: Path(p) for p in request['tracks']}
    if len(track_paths) != len(request['tracks']): raise ValueError('Duplicate track filenames')
    expected = {f'{aid}.{t}.gff3' for t in TRACKS if track_state['tracks'][t]['status'] == 'ready'}
    if set(track_paths) != expected: raise ValueError('Missing/extra selected genome tracks')
    for t in TRACKS:
        s = track_state['tracks'][t]
        if s['status'] not in ('ready', 'disabled'): raise ValueError('Unaccepted genome track')
        if s['status'] == 'ready':
            if s['reference_identity'] != ref or s['gff3_sha256'] != sha(track_paths[f'{aid}.{t}.gff3']):
                raise ValueError('Selected genome track changed: '+t)
    coding = identified/f'{aid}.gff3'; functional = request.get('functional', {})
    if request['functional_enabled']:
        fs = accepted(functional['status'], aid, stage='functional_annotation')
        if fs['original_gff3_sha256'] != sha(coding) or fs['annotated_gff3_sha256'] != sha(functional['gff']) or fs['table_sha256'] != sha(functional['table']):
            raise ValueError('Functional annotation provenance mismatch')
        if geometry(coding) != geometry(functional['gff']): raise ValueError('Functional decoration changed coding models')
        coding = Path(functional['gff'])
    elif functional: raise ValueError('Unexpected functional results when disabled')
    tx = models(coding, lengths)
    for suffix in ('proteins', 'cds'):
        if set(sequences(identified/f'{aid}_{suffix}.fa')) != set(tx): raise ValueError('Release coding FASTA IDs differ')
    def table_ids(path, column):
        with Path(path).open() as h: return [r[column] for r in csv.DictReader(h, delimiter='\t')]
    genes = {m['gene_id'] for m in tx.values()}
    assessed = table_ids(identified/'transcript_assessment.tsv', 'transcript_id')
    if len(assessed) != len(set(assessed)) or set(assessed) != set(tx): raise ValueError('Released assessment IDs differ')
    if request['functional_enabled']:
        ids = table_ids(functional['table'], 'protein_id')
        if len(ids) != len(set(ids)) or set(ids) != set(tx): raise ValueError('Functional table IDs differ from coding models')
        ids = table_ids(functional['genes'], 'gene_id')
        if len(ids) != len(set(ids)) or set(ids) != genes: raise ValueError('Gene product table IDs differ')
    # Complete upstream maps retain excluded decisions; release does not renumber.
    for name, column, ids in [('gene_id_map.tsv', 'final_gene_id', genes), ('transcript_id_map.tsv', 'final_transcript_id', set(tx))]:
        mapped = {v for v in table_ids(identified/name, column) if v not in ('', '.')}
        if mapped != ids: raise ValueError('Final ID map does not cover released models')
    inputs = {}; outputs = {}
    def copy(source, relative):
        destination = out/relative; destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        inputs[relative] = sha(source)
    copy(request['genome'], f'{aid}.reference.fa'); copy(request['masked'], f'{aid}.softmasked.fa')
    copy(request['identity'], 'provenance/reference.json')
    copy(identified/'annotation_status.json', 'provenance/annotation_status.json')
    copy(request['repeat_status'], 'provenance/reference_repeat_status.json')
    copy(request['track_status'], 'provenance/genome_tracks_status.json')
    for suffix in ('proteins', 'cds'): copy(identified/f'{aid}_{suffix}.fa', f'{aid}_{suffix}.fa')
    for name in ('gene_id_map.tsv', 'transcript_id_map.tsv', 'feature_id_map.tsv'): copy(identified/name, 'id_maps/'+name)
    copy(identified/'transcript_assessment.tsv', 'tables/transcript_assessment.tsv')
    if request['functional_enabled']:
        for key, dest in [('table', 'functional_annotations.tsv'), ('genes', 'gene_product_summary.tsv')]: copy(functional[key], 'tables/'+dest)
        copy(functional['status'], 'provenance/functional_status.json')
    owners = {}; counts = {}; selected_trnas = None
    combined = out/f'{aid}.gff3'
    write_gff(combined, [], lengths)
    map_handle = (out/'id_maps/track_id_map.tsv').open('w', newline='')
    map_writer = csv.DictWriter(map_handle, fieldnames=['track', 'original_id', 'release_id'], delimiter='\t')
    map_writer.writeheader()
    sources = [('coding', coding), ('repeats', request['repeats'])]+[(t, track_paths[f'{aid}.{t}.gff3']) for t in TRACKS if f'{aid}.{t}.gff3' in track_paths]
    try:
        for track, source in sources:
            rows = gff_records(source, lengths); inputs[f'source:{track}'] = sha(source)
            for c, a in rows:
                old = a['ID']
                if track == 'repeats':
                    a['ID'] = quote(aid+'.repeat.'+unquote(old), safe='._:-')
                    for key in ('Parent', 'Derives_from'):
                        if key in a: a[key] = ','.join(quote(aid+'.repeat.'+unquote(p), safe='._:-') for p in a[key].split(','))
                if a['ID'] in owners and owners[a['ID']] != track: raise ValueError('ID collision between release tracks: '+a['ID'])
                owners[a['ID']] = track
                map_writer.writerow(dict(track=track, original_id=old, release_id=a['ID']))
            counts[track] = dict(features=len(rows), feature_types=dict(Counter(c[2] for c, a in rows)))
            if track in TRACKS and len(rows) != track_state['tracks'][track]['features']: raise ValueError('Selected track count differs')
            if track == 'coding':
                target = out/f'{aid}.coding.gff3'
            else:
                (out/'tracks').mkdir(exist_ok=True); target = out/f'tracks/{aid}.{track}.gff3'
            write_gff(target, rows, lengths)
            if track == 'trnascan': selected_trnas = trna_models(rows)
            # Each track is already validated, with disjoint IDs across tracks.
            # Stream the validated payloads rather than retaining the entire genome twice.
            with combined.open('a') as destination, target.open() as source_handle:
                for line in source_handle:
                    if not line.startswith('#'): destination.write(line)
            del rows
            gc.collect()
    finally:
        map_handle.close()
    exports = [('transcripts', tx, out/f'{aid}_transcripts.fa')]
    if selected_trnas is not None: exports.append(('trna', selected_trnas, out/f'{aid}_tRNA.fa'))
    sequence_counts = extract_spliced(request['genome'], exports)
    lines = [f'# Annotation release: {aid}', '', 'Status: **ready**', '',
        f'Protein-coding genes: {len(genes):,}. Coding transcripts: {len(tx):,}. Functional annotation: {"included" if request["functional_enabled"] else "disabled"}.', '',
        f'- [Combined annotation]({aid}.gff3): coding models, repeats and enabled selected genome tracks.',
        f'- [Coding annotation]({aid}.coding.gff3), [proteins]({aid}_proteins.fa), [CDS]({aid}_cds.fa), [spliced transcripts]({aid}_transcripts.fa).',
        f'- [Reference]({aid}.reference.fa), [softmasked reference]({aid}.softmasked.fa).',
        '- [Complete gene ID map](id_maps/gene_id_map.tsv), [transcript ID map](id_maps/transcript_id_map.tsv), [feature ID map](id_maps/feature_id_map.tsv), [track ID map](id_maps/track_id_map.tsv).',
        '- [Transcript assessment](tables/transcript_assessment.tsv).', '']
    if selected_trnas is not None: lines += [f'[Selected spliced tRNAs]({aid}_tRNA.fa): {len(selected_trnas):,} sequences.', '']
    lines += ['Transcript sequences join the final selected exons, including accepted UTRs. tRNA sequences join annotated exons, or use the tRNA span when no exons are present. Both exports use the unmasked reference, transcript orientation and final feature IDs; headers contain identifiers only. Disabled tRNA tracks have no sequence deliverable; enabled no-hit tracks have an empty FASTA.', '']
    if request['functional_enabled']: lines += ['[Functional annotations](tables/functional_annotations.tsv) · [Gene products](tables/gene_product_summary.tsv)', '']
    lines += ['| Track | Status | Feature records | File |', '|---|---|---:|---|']
    for t in ('repeats',)+TRACKS:
        enabled = t in counts
        lines.append(f'| {t} | {"ready" if enabled else "disabled"} | {counts[t]["features"] if enabled else "—"} | '+(f'[{t}](tracks/{aid}.{t}.gff3)' if enabled else '—')+' |')
    lines += ['', 'Counts include child and multipart records. Overlapping tracks remain separate evidence; they are not deduplicated genes. Coding IDs, sequences and geometry are unchanged. Repeat IDs are namespaced to this assembly; the track map records the conversion.', '',
        'The combined file is grouped by track with parents before children, not globally coordinate sorted. External evidence IDs in Target/other annotations are preserved.', '',
        '[Release manifest](release_status.json) · [Checksums](checksums.sha256) · [Source manifests](provenance/)', '',
        'Detailed run reports and raw/audit results (links require the surrounding run directory):', '',
        f'- [Structural QC](../../qc/structural/{aid}/README.md)',
        f'- [Genome tracks, raw calls and filter decisions](../../genome_tracks/{aid}/README.md)',
        f'- [Model selection](../../structural_annotation/selection/{aid}/README.md)',
        '- [Reference and repeat preparation](../../reference_repeats/README.md)', '']
    if request['functional_enabled']: lines += [f'- [Functional searches and provenance](../../functional_annotation/{aid}/README.md)', '']
    (out/'README.md').write_text('\n'.join(lines), encoding='utf-8')
    for path in sorted(out.rglob('*')):
        if path.is_file(): outputs[path.relative_to(out).as_posix()] = dict(sha256=sha(path), bytes=path.stat().st_size)
    state = dict(schema_version=1, stage='annotation_release', status='ready', assembly_id=aid, reference_identity=ref,
        genes=len(genes), coding_transcripts=len(tx), sequence_counts=sequence_counts, functional_enabled=request['functional_enabled'], tracks=counts,
        disabled_tracks=[t for t in TRACKS if t not in counts], source_sha256=inputs, outputs=outputs,
        scope='Validated final deliverables; no model changes, new filtering or searches')
    dump(out/'release_status.json', state)
    with (out/'checksums.sha256').open('w') as h:
        for name, value in sorted(outputs.items()): h.write(value['sha256']+'  '+name+'\n')
        h.write(sha(out/'release_status.json')+'  release_status.json\n')
    return state


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True); parser.add_argument('--outdir', required=True)
    args = parser.parse_args(); package(json.loads(base64.b64decode(args.request)), args.outdir)
