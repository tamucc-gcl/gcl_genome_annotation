#!/usr/bin/env python3
"""Alignment reporting only; report changes do not invalidate mapping."""
import argparse
import subprocess
from pathlib import Path
from rna_alignment import decode, dump, load, report


def alignment_multiqc(request, outdir):
    """Summarize native staged logs, keeping assembly/readset identities distinct."""
    folder = Path(outdir) / 'multiqc'
    folder.mkdir(parents=True, exist_ok=True)
    paths, samples = request['paths'], request['samples']
    if len(paths) != len(samples):
        raise ValueError('Alignment MultiQC log/sample count mismatch')
    if not paths:
        (folder / 'multiqc_report.html').write_text(
            '<!doctype html><html><head><title>RNA alignment QC</title></head>'
            '<body><h1>RNA alignment QC</h1><p>No RNA alignment was requested.</p></body></html>\n',
            encoding='utf-8')
        return
    names = {Path(path).stem: sample for path, sample in zip(paths, samples)}
    if len(names) != len(paths):
        raise ValueError('Alignment MultiQC staged log names must be unique')
    config = folder / 'multiqc_config.json'
    dump(config, {'title': 'RNA alignment QC', 'use_filename_as_sample_name': True,
                  'fn_clean_exts': ['.txt'], 'sample_names_replace_exact': True,
                  'sample_names_replace': names})
    inputs = folder / 'input_files.txt'
    inputs.write_text(''.join(str(Path(path).resolve()) + '\n' for path in paths), encoding='utf-8')
    with (folder / 'multiqc.log').open('w', encoding='utf-8') as log:
        subprocess.run(['multiqc', '--file-list', str(inputs), '--config', str(config),
                        '--module', 'bowtie2', '--module', 'samtools',
                        '--require-logs', '--force',
                        '--outdir', str(folder), '--filename', 'multiqc_report.html'],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    if not (folder / 'multiqc_report.html').is_file():
        raise RuntimeError('MultiQC produced no RNA alignment report; inspect multiqc.log')



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'qc', 'strands', 'request', 'multiqc-request', 'outdir'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    alignment_multiqc(decode(args.multiqc_request), args.outdir)
    result = report(load(args.plan), [load(x) for x in Path(args.qc).glob('*.json')],
                    [load(x) for x in Path(args.strands).glob('*.json')],
                    decode(args.request))
    dump(Path(args.outdir) / 'rna_alignment_status.json', result)
    lines = ['# RNA alignment', '', f'Status: **{result["status"]}**', '',
             'Mapping uses the accepted normalized reference and previously selected reads.',
             'Automatic strandedness uses unique, high-quality spliced primary alignments (R1 for paired reads).',
             'Validated readset BAMs are listed for BRAKER; gene prediction has not run.', '',
             '| Assembly | Readset | Mapped | Informative | Strand | Status |', '|---|---|---:|---:|---|---|']
    for row in result['readsets']:
        state = row['strand']
        lines.append(f'| {state["assembly_id"]} | {state["readset_id"]} | {row.get("mapped_fraction", 0):.1%} | {state["informative"]} | {state.get("strandedness", "unresolved")} | {state["status"]} |')
        if state.get('reason'):
            lines.extend(['', state['reason'], ''])
    lines.extend(['', '| Assembly | Library | Status | Strand | Detail |', '|---|---|---|---|---|'])
    for row in result['libraries']:
        detail = ', '.join(f'[BAM {i}]({path})' for i, path in enumerate(row['bams'], 1)) if row['status'] == 'ready' else row['reason']
        lines.append(f'| {row["assembly_id"]} | {row["library_id"]} | {row["status"]} | {row.get("strandedness", "unresolved")} | {detail} |')
    lines.extend(['', '[Alignment QC: HISAT2 and samtools](multiqc/multiqc_report.html)',
                  '[Resolved alignment statistics and BAM inputs](rna_alignment_status.json)', ''])
    Path(args.outdir, 'README.md').write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
