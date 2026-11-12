#!/usr/bin/env python3
"""Export centroAnno satellite-array candidates using verified coordinate conventions."""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from reference_repeats import sha as checksum

def monomers(path, convention, lengths):
    """Require explicit schema/convention. New C++ output is 0-based inclusive."""
    if convention not in ('zero_inclusive','one_inclusive'):raise ValueError('Explicit coordinate convention required')
    result=[]
    with Path(path).open() as h:
        for c in csv.reader(h):
            if not c or c[0].strip().lower()=='sequence name':continue
            if len(c)==7:seq,name,strand,begin,end,identity,length=c
            elif len(c)==6:
                seq,name,begin,end,identity,length=c;strand='-' if name.endswith("'") else '+';name=name.rstrip("'")
            else:raise ValueError('Unknown monomer schema')
            a,b=int(begin),int(end);n=int(length)
            if b-a+1!=n:raise ValueError('Monomer coordinate/length mismatch')
            if convention=='one_inclusive':a-=1;b-=1
            if seq not in lengths or not 0<=a<=b<lengths[seq]:raise ValueError('Monomer outside reference')
            if strand not in ('+','-') or not 0<=float(identity)<=1:raise ValueError('Invalid monomer strand/identity')
            result.append(dict(seq=seq,name=name,strand=strand,start=a,end=b+1,identity=float(identity),length=n))
    return result


def patch_centrofinder(source, output):
    """Correction of inclusive monomer ends entering half-open arrays."""
    text=Path(source).read_text()
    replacements={
        "m['end'] > c['start']":"m['end'] + 1 > c['start']",
        "block_bp += region_monos[j]['end'] - region_monos[j]['start']":"block_bp += region_monos[j]['end'] - region_monos[j]['start'] + 1",
        "new_end = region_monos[max(j, 0)]['end']":"new_end = region_monos[max(j, 0)]['end'] + 1",
        "new_end = region_monos[i]['end']":"new_end = region_monos[i]['end'] + 1"}
    for old,new in replacements.items():
        expected=2 if old.startswith('block_bp') else 1
        if text.count(old)!=expected:raise ValueError('Pinned centroFinder patch context changed: '+old)
        text=text.replace(old,new)
    with Path(output).open('x') as h:h.write(text)
    return {'native_sha256':checksum(Path(source)),'experimental_sha256':checksum(Path(output)),
            'changes':replacements,'purpose':'Inclusive monomer end to half-open candidate-array boundary'}



def run(aid, genome, identity, tool, threads, min_monomer, gff, outdir):
    out = Path(outdir)
    out.mkdir()
    tool = Path(tool).resolve()
    lengths = {r['id']:r['length'] for r in json.loads(Path(identity).read_text())['contigs']}
    patched = out/'centroFinder_halfopen.py'
    patch = patch_centrofinder(tool/'script/centroFinder.py', patched)
    seen = {}
    with Path(gff).open('x') as gff_handle, (out/'monomers.tsv').open('x',newline='') as monomers_handle, \
         (out/'hors.tsv').open('x',newline='') as hors_handle, (out/'contigs.tsv').open('x',newline='') as contigs_handle:
        gff_handle.write('##gff-version 3\n')
        monomer_writer = csv.DictWriter(monomers_handle, ['seq','name','strand','start','end','identity','length'],delimiter='\t')
        hor_writer = csv.writer(hors_handle,delimiter='\t')
        contig_writer = csv.writer(contigs_handle,delimiter='\t')
        monomer_writer.writeheader()
        hor_writer.writerow(['contig','start','end','native'])
        contig_writer.writerow(['contig','length','monomer_count','array_count','localization_status'])
        for ordinal, record in enumerate(SeqIO.parse(str(genome),'fasta'),1):
            if record.id in seen or lengths.get(record.id) != len(record):
                raise ValueError('Reference sequence differs from contig dictionary')
            seen[record.id] = len(record)
            folder = out/f'contig_{ordinal:04d}'
            folder.mkdir()
            fasta = folder/'input.fa'
            SeqIO.write(record,str(fasta),'fasta')
            raw = folder/'native'
            raw.mkdir()
            subprocess.run([tool/'centroAnno',fasta,'-o',raw,'-x','anno-asm','-t',str(threads)],check=True)
            csvs = list(raw.rglob('*_decomposedResult.csv'))
            if len(csvs) > 1:
                raise ValueError('Ambiguous decomposition per contig')
            calls = monomers(csvs[0],'zero_inclusive',{record.id:len(record)}) if csvs else []
            fasta.unlink()
            count = 0
            status = 'no_monomer_calls'
            if calls:
                monomer_writer.writerows(calls)
                prefix = folder/'hor'
                subprocess.run([sys.executable,tool/'script/detectHORs.py','-i',csvs[0],'-o',prefix],check=True)
                with Path(str(prefix)+'_allHORs.csv').open() as handle:
                    for row in csv.DictReader(handle):
                        start,end = int(row['Start position']),int(row['End position'])+1
                        if row['Sequence name'] != record.id or not 0 <= start < end <= len(record):
                            raise ValueError('HOR outside reference')
                        hor_writer.writerow([record.id,start,end,json.dumps(row,sort_keys=True)])
                status = 'below_localization_window'
                if len(record) >= 50000:
                    SeqIO.write(record,str(fasta),'fasta')
                    prefix = folder/'arrays'
                    subprocess.run([sys.executable,patched,'--mono-csv',csvs[0],'--fasta',fasta,'--out-prefix',prefix,
                        '--no-offset','--min-mono-len',str(min_monomer)],check=True)
                    fasta.unlink()
                    status = 'assessed'
                    with Path(str(prefix)+'_censatarray_candidates.bed').open() as handle:
                        for line in handle:
                            if not line.strip() or line.startswith('#'):
                                continue
                            seq,start,end,*rest = line.rstrip().split('\t')
                            start,end = int(start),int(end)
                            if seq != record.id or not 0 <= start < end <= len(record):
                                raise ValueError('Candidate array outside reference')
                            count += 1
                            attrs = f'ID=candidate_array.{ordinal:04d}.{count:03d};Name=candidate_satellite_array;functional_centromere=unvalidated;native_rank={count}'
                            gff_handle.write(f'{seq}\tcentroAnno_candidate\trepeat_region\t{start+1}\t{end}\t.\t.\t.\t{attrs}\n')
            contig_writer.writerow([record.id,len(record),len(calls),count,status])
    if seen != lengths:
        raise ValueError('Reference is missing contigs')
    (out/'settings.json').write_text(json.dumps({'assembly_id':aid,'minimum_monomer_length':min_monomer,
        'centroanno_revision':(tool/'revision.txt').read_text().strip(),
        'table_coordinates':'0-based half-open','patch':patch},indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('aid','genome','identity','tool','gff','outdir'):
        parser.add_argument('--'+name,required=True)
    parser.add_argument('--threads',type=int,required=True)
    parser.add_argument('--min-monomer',type=int,required=True)
    args = parser.parse_args()
    if args.threads < 1 or args.min_monomer < 1:
        parser.error('threads and minimum monomer length must be positive')
    run(args.aid,args.genome,args.identity,args.tool,args.threads,args.min_monomer,args.gff,args.outdir)
