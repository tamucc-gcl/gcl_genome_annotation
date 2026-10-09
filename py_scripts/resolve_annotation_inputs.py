#!/usr/bin/env python3
"""Offline annotation preflight. No acquisition or annotation tools are launched."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ANALYSES = ('busco', 'omark', 'eggnog', 'diamond', 'interproscan', 'kofamscan',
            'deepsig', 'infernal', 'trnascan', 'telomere', 'centromere')
SPECIES_FIELDS = ('species', 'orthodb_clade', 'busco_lineage', 'taxa_name',
                  'telomere_motif', 'genetic_code', 'eggnog_taxid')
ASSEMBLY_FIELDS = set(('assembly_id', 'assembly_fasta', 'assembly_path', 'taxid',
    'ncbi_speciesID', 'individual_id', 'haplotype', 'comparison_group', 'rna_readsets',
    'use_protein_evidence', 'protein_fastas', 'ref_proteome_accession',
    'min_contig_length') + SPECIES_FIELDS + tuple('run_' + x for x in ANALYSES))
RNA_FIELDS = set(('readset_id', 'library_id', 'sample_id', 'tissue', 'layout',
    'fastq_1', 'fastq_2', 'run_accession', 'strandedness', 'trim_adapters',
    'trim_quality', 'subset_fraction', 'subset_seed'))
DEFAULTS = {'use_protein_evidence': True, 'min_contig_length': 1000,
            'trim_adapters': True, 'trim_quality': True, 'subset_fraction': 1.0,
            'subset_seed': 42, **{'run_' + x: True for x in ANALYSES}}
ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
RUN = re.compile(r'(SRR|ERR|DRR)[0-9]+\Z')


def present(value):
    return value is not None and str(value).strip() != ''


def selected(row, params, key, default=None):
    """Do not collapse explicit False, 0 or 'none' into an inherited value."""
    if present(row.get(key)):
        return row[key], 'assembly_row'
    if present(params.get(key)):
        return params[key], 'global_parameter'
    return default, 'default'


def boolean(value):
    token = str(value).strip().lower()
    if token not in ('true', 'false'):
        raise ValueError(f'Expected true/false, received {value!r}')
    return token == 'true'


def integer(value, minimum=0):
    if not re.fullmatch(r'[0-9]+', str(value).strip()):
        raise ValueError(f'Expected integer >= {minimum}, received {value!r}')
    result = int(value)
    if result < minimum:
        raise ValueError(f'Expected integer >= {minimum}, received {value!r}')
    return result


def identifier(value):
    if not ID.fullmatch(value or '') or value in ('.', '..'):
        raise ValueError(f'Invalid identifier {value!r}')
    return value


def items(value):
    if not present(value) or str(value).lower() == 'none':
        return []
    result = [x.strip() for x in str(value).split(';')]
    if any(not x for x in result) or len(set(result)) != len(result):
        raise ValueError('Empty or duplicate list item')
    return result


def path_at(value, base):
    candidate = Path(value).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (base / candidate).resolve()


def readable(value, base):
    path = path_at(value, base)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty file: {path}')
    with path.open('rb') as handle:
        handle.read(1)
    return str(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_table(path, allowed, required):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        header = [x.strip() for x in header]
        if not header or any(not x for x in header) or len(set(header)) != len(header):
            raise ValueError(f'{path}: missing/duplicate/empty header')
        if not set(required) <= set(header):
            raise ValueError(f'{path}: required columns: {required}')
        unknown = set(header) - allowed
        if unknown:
            raise ValueError(f'{path}: unknown columns: {sorted(unknown)}')
        rows = []
        for cells in reader:
            if not any(x.strip() for x in cells):
                continue
            if len(cells) != len(header):
                raise ValueError(f'{path}:{reader.line_num}: expected {len(header)} fields')
            rows.append(dict(zip(header, (x.strip() for x in cells))))
    return rows


def alias(row, canonical, old):
    row = dict(row)
    if present(row.get(old)):
        if present(row.get(canonical)) and str(row[canonical]) != str(row[old]):
            raise ValueError(f'Conflicting {canonical} and legacy alias {old}')
        row[canonical] = row[old]
    return row


class Taxonomy:
    """Read a fixed local NCBI taxdump; ancestry IDs match a versioned rules file."""
    def __init__(self, directory):
        self.nodes, self.names, self.merged = {}, {}, {}
        self.checksums = {}
        directory = Path(directory)
        for filename in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
            path = directory / filename
            if filename == 'merged.dmp' and not path.exists():
                continue
            self.checksums[filename] = digest(path)
            with path.open(encoding='utf-8') as handle:
                for line in handle:
                    fields = [x.strip() for x in line.split('|')]
                    needed = {'nodes.dmp': 8, 'names.dmp': 4, 'merged.dmp': 2}[filename]
                    if len(fields) < needed:
                        raise ValueError(f'Malformed taxonomy record in {filename}')
                    if filename == 'nodes.dmp':
                        self.nodes[fields[0]] = (fields[1], int(fields[6]), fields[7] == '1')
                    elif filename == 'names.dmp' and fields[3] == 'scientific name':
                        self.names[fields[0]] = fields[1]
                    elif filename == 'merged.dmp':
                        self.merged[fields[0]] = fields[1]
        self.cache = {}

    def resolve(self, requested):
        taxid = str(integer(requested, 1))
        if taxid in self.cache:
            return self.cache[taxid]
        current, seen = taxid, set()
        while current in self.merged:
            if current in seen:
                raise ValueError(f'Cycle resolving merged taxid {taxid}')
            seen.add(current)
            current = self.merged[current]
        canonical, ancestors, seen, code = current, [], set(), None
        while True:
            if current in seen or current not in self.nodes:
                raise ValueError(f'Unknown/broken taxonomy lineage for taxid {taxid}')
            seen.add(current)
            ancestors.append(current)
            parent, node_code, inherited = self.nodes[current]
            if code is None and not inherited:
                code = node_code
            if current == parent:
                break
            current = parent
        if '2759' not in ancestors:
            raise ValueError(f'taxid {taxid} is outside the supported eukaryotic scope')
        if canonical not in self.names or not code:
            raise ValueError(f'Incomplete taxonomy record for taxid {taxid}')
        result = {'requested_taxid': taxid, 'taxid': canonical,
                  'species': self.names[canonical], 'ancestors': ancestors,
                  'lineage': [self.names.get(x, x) for x in reversed(ancestors)],
                  'genetic_code': code}
        self.cache[taxid] = result
        return result


def derive(key, taxon, rules):
    if not taxon:
        return None
    if key == 'species':
        return taxon['species']
    if key == 'genetic_code':
        return taxon['genetic_code']
    if key == 'taxa_name':
        return taxon['taxid']  # Dfam accepts a taxonomy ID; no invented clade name.
    # eggNOG taxonomic restriction is a deliberate biological override, not species identity.
    if key == 'eggnog_taxid':
        return 'none'
    options = {str(x['taxid']): x['value'] for x in rules.get(key, [])}
    return next((options[x] for x in taxon['ancestors'] if x in options), None)


def resolve_assembly(row, params, tax, rules, sheet_base, launch_base):
    row = alias(alias(row, 'taxid', 'ncbi_speciesID'), 'assembly_fasta', 'assembly_path')
    result = {'assembly_id': identifier(row['assembly_id']), 'settings': {}, 'origins': {}}
    result['assembly_fasta'] = readable(row.get('assembly_fasta', ''), sheet_base)
    result['rna_readsets'] = sorted(items(row.get('rna_readsets')))
    for key in ('individual_id', 'haplotype', 'comparison_group'):
        result[key] = row.get(key) or None
    settings, origins = result['settings'], result['origins']
    for key in ('use_protein_evidence', 'min_contig_length') + tuple('run_' + x for x in ANALYSES):
        value, origins[key] = selected(row, params, key, DEFAULTS[key])
        settings[key] = integer(value) if key == 'min_contig_length' else boolean(value)
    taxid, origins['taxid'] = selected(row, params, 'taxid')
    taxon = None
    if present(taxid):
        if tax is None:
            raise ValueError('taxid supplied but no local --taxdump_dir was provided')
        taxon = tax.resolve(taxid)
    settings['taxid'] = taxon['taxid'] if taxon else None
    result['taxonomy'] = taxon
    for key in SPECIES_FIELDS:
        value, origin = selected(row, params, key, 'auto')
        if str(value).lower() == 'auto':
            value, origin = derive(key, taxon, rules), 'taxonomy' if taxon else 'unresolved'
            if key == 'eggnog_taxid':
                value, origin = 'none', 'default_unrestricted'
        settings[key], origins[key] = value, origin
    if not present(settings['genetic_code']):
        raise ValueError('Cannot derive nuclear genetic_code; provide taxid/taxdump or an explicit override')
    settings['genetic_code'] = integer(settings['genetic_code'], 1)
    for key in ('protein_fastas', 'ref_proteome_accession'):
        value, origin = selected(row, params, key)
        sources = items(value)
        if not settings['use_protein_evidence']:
            if present(row.get(key)) and sources:
                raise ValueError(f'{key} explicitly supplied while protein evidence is disabled')
            sources = []
        if key == 'protein_fastas':
            base = sheet_base if origin == 'assembly_row' else launch_base
            sources = [readable(x, base) for x in sources]
        else:
            for accession in sources:
                if not re.fullmatch(r'GC[AF]_[0-9]+\.[0-9]+', accession):
                    raise ValueError(f'Expected versioned NCBI assembly accession: {accession}')
        settings[key], origins[key] = sources, origin
    clade = settings['orthodb_clade']
    if not settings['use_protein_evidence']:
        if present(row.get('orthodb_clade')) and row['orthodb_clade'].lower() != 'none':
            raise ValueError('Explicit OrthoDB request conflicts with disabled protein evidence')
        settings['orthodb_clade'] = 'none'
    elif not clade:
        raise ValueError('Cannot derive OrthoDB partition: correct taxid or override orthodb_clade')
    elif clade != 'none' and clade not in rules['allowed_orthodb_clades']:
        raise ValueError(f'Unsupported OrthoDB partition {clade!r}')
    for field, enabled in (('busco_lineage', settings['run_busco']),
                           ('telomere_motif', settings['run_telomere']),
                           ('taxa_name', True)):
        if enabled and (not present(settings[field]) or settings[field] == 'none'):
            raise ValueError(f'Cannot derive {field}; supply an override or disable its optional analysis')
    if settings['run_omark'] and not settings['taxid']:
        raise ValueError('OMArk requires an effective taxid; supply taxid or disable run_omark')
    if settings['telomere_motif'] not in (None, 'none'):
        motif = str(settings['telomere_motif']).upper()
        if not re.fullmatch('[ACGT]+', motif):
            raise ValueError('telomere_motif must contain only A/C/G/T')
        settings['telomere_motif'] = motif
    proteins = settings['use_protein_evidence'] and bool(
        settings['orthodb_clade'] != 'none' or settings['protein_fastas'] or settings['ref_proteome_accession'])
    if settings['use_protein_evidence'] and not proteins:
        raise ValueError('Protein evidence enabled but every source is disabled')
    result['prediction_mode'] = ('rna_protein' if proteins else 'rna') if result['rna_readsets'] else (
        'protein' if proteins else 'ab_initio')
    result['status'] = 'ready'
    return result


def resolve_rna(row, params, base):
    result = {k: row.get(k) or None for k in RNA_FIELDS}
    for key in ('readset_id', 'library_id'):
        identifier(row.get(key))
    if row.get('layout') not in ('paired', 'single'):
        raise ValueError('layout must be paired or single')
    accession = row.get('run_accession')
    if accession:
        if not RUN.fullmatch(accession) or row.get('fastq_1') or row.get('fastq_2'):
            raise ValueError('Supply one SRR/ERR/DRR run accession OR local FASTQ, not both')
    else:
        if not row.get('fastq_1') or bool(row.get('fastq_2')) != (row['layout'] == 'paired'):
            raise ValueError('Read files do not match declared layout')
        for key in ('fastq_1', 'fastq_2'):
            if row.get(key):
                result[key] = readable(row[key], base)
        if result['fastq_1'] == result['fastq_2']:
            raise ValueError('Paired mates refer to the same file')
    strand, _ = selected(row, params, 'strandedness', 'auto')
    if strand not in ('auto', 'unstranded', 'forward', 'reverse'):
        raise ValueError('Invalid strandedness')
    result['strandedness'] = strand
    for key in ('trim_adapters', 'trim_quality'):
        result[key] = boolean(selected(row, params, key, DEFAULTS[key])[0])
    result['subset_fraction'] = float(selected(row, params, 'subset_fraction', 1.0)[0])
    if not 0 < result['subset_fraction'] <= 1:
        raise ValueError('subset_fraction must be >0 and <=1')
    result['subset_seed'] = integer(selected(row, params, 'subset_seed', 42)[0])
    result['status'] = 'ready'
    result['content_validation'] = 'pending_cluster_read_QC'
    return result


def build_manifest(request, rules):
    params = alias(request.get('params', {}), 'taxid', 'ncbi_speciesID')
    launch = Path(request['launch_dir']).resolve()
    assembly_sheet = path_at(request['sample_sheet'], launch)
    rows = read_table(assembly_sheet, ASSEMBLY_FIELDS, ['assembly_id'])
    if not rows:
        raise ValueError('Assembly sheet has no rows')
    rna_sheet = path_at(request['rna_sheet'], launch) if request.get('rna_sheet') else None
    rna_rows = read_table(rna_sheet, RNA_FIELDS, ['readset_id', 'library_id', 'layout']) if rna_sheet else []
    for table, key in ((rows, 'assembly_id'), (rna_rows, 'readset_id')):
        ids = [identifier(row.get(key)) for row in table]
        if len(ids) != len(set(ids)):
            raise ValueError(f'Duplicate {key}')
    taxdir = request.get('taxdump_dir')
    tax = Taxonomy(path_at(taxdir, launch)) if taxdir else None
    rna, owners, libraries = {}, {}, {}
    for row in rna_rows:
        rid = row['readset_id']
        try:
            rna[rid] = resolve_rna(row, params, rna_sheet.parent)
        except (ValueError, OSError) as exc:
            rna[rid] = {'readset_id': rid, 'status': 'failed', 'reason': str(exc)}
            continue
        resolved = rna[rid]
        for source in ('fastq_1', 'fastq_2', 'run_accession'):
            value = resolved.get(source)
            if value and value in owners:
                raise ValueError(f'RNA source assigned twice: {value}; reuse readset {owners[value]}')
            if value:
                owners[value] = rid
        signature = tuple(resolved[k] for k in ('layout', 'strandedness', 'trim_adapters', 'trim_quality'))
        lib = resolved['library_id']
        if lib in libraries and libraries[lib] != signature:
            raise ValueError(f'Inconsistent settings within RNA library {lib}')
        libraries[lib] = signature
    assemblies = []
    for row in rows:
        try:
            resolved = resolve_assembly(row, params, tax, rules, assembly_sheet.parent, launch)
            for rid in resolved['rna_readsets']:
                if rid not in rna or rna[rid]['status'] != 'ready':
                    raise ValueError(f'Required RNA read set {rid!r} is unknown or failed validation')
        except (ValueError, OSError) as exc:
            resolved = {'assembly_id': row['assembly_id'], 'status': 'blocked', 'reason': str(exc)}
        assemblies.append(resolved)
    return {'schema_version': 1, 'stage': 'input_preflight',
            'validation_scope': 'settings_and_file_presence_only',
            'pending_checks': ['read_content_and_accession_acquisition', 'tool_dataset_compatibility',
                               'reference_sequence_integrity', 'strandedness_inference'],
            'status': 'ready' if all(x['status'] == 'ready' for x in assemblies) and all(
                x['status'] == 'ready' for x in rna.values()) else 'incomplete',
            'assemblies': assemblies, 'rna_readsets': list(rna.values()),
            'taxonomy_checksums': tax.checksums if tax else {},
            'input_sheet_checksums': {'assembly': digest(assembly_sheet),
                                     'rna': digest(rna_sheet) if rna_sheet else None}}


def write_report(manifest, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'resolved_inputs.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    def escape(value):
        return str(value).replace('|', '\\|').replace('\n', ' ').replace('\r', ' ')
    lines = ['# Annotation input preflight', '', f"Status: **{manifest['status']}**", '',
             'This report validates configuration and file presence, not read contents or biological results.', '',
             '| Assembly | Status | Mode | Reason |', '|---|---|---|---|']
    for row in manifest.get('assemblies', []):
        lines.append('| ' + ' | '.join(escape(row.get(k, '')) for k in (
            'assembly_id', 'status', 'prediction_mode', 'reason')) + ' |')
    for row in manifest.get('rna_readsets', []):
        if row['status'] != 'ready':
            lines.extend(['', f"RNA `{escape(row['readset_id'])}`: {escape(row['reason'])}"])
    if manifest.get('error'):
        lines.extend(['', escape(manifest['error'])])
    lines.extend(['', '[Resolved settings and provenance](resolved_inputs.json)', ''])
    (outdir / 'README.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True)
    parser.add_argument('--rules', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--report-only', action='store_true', help='Write failures as data for Nextflow status assertion')
    args = parser.parse_args()
    try:
        request = json.loads(Path(args.request).read_text(encoding='utf-8'))
        rules = json.loads(Path(args.rules).read_text(encoding='utf-8'))
        manifest = build_manifest(request, rules)
        manifest['rules_sha256'] = digest(args.rules)
    except (ValueError, OSError, KeyError, csv.Error) as exc:
        manifest = {'schema_version': 1, 'status': 'failed', 'error': str(exc)}
    write_report(manifest, Path(args.outdir))
    return 0 if args.report_only or manifest['status'] == 'ready' else 2


if __name__ == '__main__':
    raise SystemExit(main())
