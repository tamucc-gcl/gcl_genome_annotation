#!/usr/bin/env python3
"""Validate Portcullis inputs and preserve native junction results, without rescoring."""
import argparse
import base64
import csv
import hashlib
import json
from pathlib import Path
import re

from prepare_protein_evidence import checksum
from reference_repeats import fasta


def decode(value):
    return json.loads(base64.b64decode(value))


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def check(genome, identity, bam, request, output):
    import pysam
    request=decode(request)
    row=request['readset']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',row['readset_id']):
        raise ValueError('Unsafe/missing readset identity')
    if row['strandedness'] not in ('forward','reverse','unstranded') or row['layout'] not in ('paired','single'):
        raise ValueError('Portcullis requires resolved library strand and layout')
    reference=json.loads(Path(identity).read_text())
    contigs=[dict(id=h.split()[0],length=len(s),uppercase_sha256=hashlib.sha256(s.upper().encode('ascii')).hexdigest())
             for h,s in fasta(genome)]
    if contigs!=reference['contigs']:
        raise ValueError('Portcullis genome differs from accepted reference')
    with pysam.AlignmentFile(str(bam),'rb') as handle:
        if dict(zip(handle.references,handle.lengths))!={r['id']:r['length'] for r in contigs}:
            raise ValueError('Portcullis BAM dictionary differs from reference')
        if handle.header.to_dict().get('HD',{}).get('SO')!='coordinate':
            raise ValueError('Portcullis requires coordinate-sorted BAMs')
        handle.check_index()
    dump(output,dict(status='ready',assembly_id=request['assembly_id'],readset=row,
        reference_identity=reference['reference_identity'],bam_bytes=Path(bam).stat().st_size,
        input_validation_scope='Accepted genome sequence identity and BAM dictionary, sort header and readable index; Nextflow manages BAM file identity/staging'))


def table(path, lengths=None):
    """Portcullis TAB start/end are 0-based inclusive; use half-open internally."""
    result={}
    with Path(path).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        required={'refname','start','end','consensus-strand','nb_raw_aln'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f'Unrecognized Portcullis table header: {path}')
        for row in reader:
            seq,start,end,strand=row['refname'],int(row['start']),int(row['end'])+1,row['consensus-strand']
            if strand not in ('+','-','?','.') or not 0<=start<end:
                raise ValueError(f'Invalid Portcullis junction: {path}')
            if lengths is not None and (seq not in lengths or end>lengths[seq]):
                raise ValueError('Portcullis junction outside accepted reference')
            count=int(row['nb_raw_aln'])
            if count<1:raise ValueError('Observed junction has no supporting alignments')
            key=(seq,start,end,strand if strand in ('+','-') else '.')
            if key in result:raise ValueError('Duplicate native junction key')
            result[key]=dict(raw_alignments=count)
    return result


def finish(folder, threshold):
    folder=Path(folder)
    state=json.loads((folder/'inputs.json').read_text())
    raw=table(folder/'all.junctions.tab')
    passed=table(folder/'filtered.pass.junctions.tab')
    failed=table(folder/'filtered.fail.junctions.tab')
    if set(passed)&set(failed) or set(raw)!=set(passed)|set(failed):
        raise ValueError('Portcullis pass/fail tables do not partition observed junctions')
    for key,value in (passed | failed).items():
        if value!=raw[key]:raise ValueError('Portcullis filtered junction counts differ from raw')
    log=(folder/'filter.log').read_text(errors='replace')
    state.update(tool='Portcullis',version=(folder/'version.txt').read_text().strip(),threshold=threshold,
        observed_junctions=len(raw),passing_junctions=len(passed),failed_junctions=len(failed),
        native_low_junction_fallback='lenient rule-based filter' in log,
        table_sha256={name:checksum(folder/name) for name in ('all.junctions.tab','filtered.pass.junctions.tab','filtered.fail.junctions.tab')})
    dump(folder/'junctions.json',state)
    return state


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('check')
    for name in ('genome','identity','bam','request','out'):p.add_argument('--'+name,required=True)
    p=sub.add_parser('finish');p.add_argument('--folder',required=True);p.add_argument('--threshold',type=float,required=True)
    args=parser.parse_args()
    if args.action=='check':check(args.genome,args.identity,args.bam,args.request,args.out)
    else:finish(args.folder,args.threshold)
