"""Small model/native-score fixtures shared by active regression tests."""
import json
from ingenannot_inputs import check
from prepare_protein_evidence import checksum

IDENTITY = {'reference_identity': 'test-ref', 'contigs': [{'id': 'chr1', 'length': 1000}]}
MODEL = '''##gff-version 3
chr1\tBRAKER\tgene\t101\t400\t.\t+\t.\tID=g1
chr1\tBRAKER\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1
chr1\tBRAKER\texon\t101\t200\t.\t+\t.\tParent=t1
chr1\tBRAKER\texon\t301\t400\t.\t+\t.\tParent=t1
chr1\tBRAKER\tCDS\t151\t200\t.\t+\t0\tParent=t1
chr1\tBRAKER\tCDS\t301\t350\t.\t+\t1\tParent=t1
'''

def native_aed_fixture(root, rna=True, protein=True):
    gff=root/'models.gff3'; gff.write_text(MODEL)
    mode='rna_protein' if rna and protein else 'rna' if rna else 'protein' if protein else 'ab_initio'
    validation=root/'validation.json'
    validation.write_text(json.dumps(dict(status='ready',stage='final_annotation',mode=mode,
        gff3_sha256=checksum(gff),reference_identity='test-ref')))
    r=root/'rna.gtf';r.write_text('# RNA fixture\n')
    p=root/'protein.gff';p.write_text('# protein fixture\n')
    identity=root/'identity.json';identity.write_text(json.dumps(IDENTITY))
    annotated=root/'native.gff3'
    annotated.write_text(MODEL.replace('ID=t1;Parent=g1',
        'ID=t1;Parent=g1;aed_ev_tr=0.3000;aed_ev_pr=0.1000;ev_tr=r1;ev_pr=p1;ev_tr_penalty=no'))
    check(gff,validation,r if rna else None,p if protein else None,root/'inputs.json')
    return gff,validation,identity,annotated
