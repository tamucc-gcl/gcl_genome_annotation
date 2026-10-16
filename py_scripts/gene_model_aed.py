"""Optional source-specific AED filtering; never interpret missing evidence as 1."""
import csv
import math


def gene_scores(path,source,maximum,genes):
    if maximum is None:
        return {}
    if source not in ('rna','protein'):
        raise ValueError('Set filter_aed_source to rna or protein when filter_max_aed is enabled')
    if not math.isfinite(maximum) or not 0<=maximum<=1:
        raise ValueError('filter_max_aed must be between 0 and 1')
    selected=source+'_aed'; result={}; seen=set()
    with open(path) as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if not {'gene_id','transcript_id',selected} <= set(reader.fieldnames or []):
            raise ValueError('AED table is missing the requested source column: '+selected)
        for row in reader:
            gid,tid=row['gene_id'],row['transcript_id']
            if gid not in genes or tid in seen:
                raise ValueError('AED table has unknown genes or duplicate transcript IDs')
            seen.add(tid)
            if row[selected] in ('','.','None'):
                raise ValueError('Requested AED source is unavailable; disable AED filtering or choose an available source')
            value=float(row[selected])
            if not math.isfinite(value) or not 0<=value<=1:
                raise ValueError('Invalid AED value')
            result[gid]=min(result.get(gid,1),value)
    if set(result)!=set(genes):
        raise ValueError('AED table does not cover all model genes')
    return result
