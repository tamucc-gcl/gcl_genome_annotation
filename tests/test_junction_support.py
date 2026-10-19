import base64
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from junction_support import report
from portcullis_inputs import check, finish, table
from prepare_protein_evidence import checksum


def encoded(value):return base64.b64encode(json.dumps(value).encode()).decode()


class JunctionSupportContracts(unittest.TestCase):
    def test_input_check_blocks_wrong_bam_dictionary_and_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            import hashlib
            genome=root/'genome.fa';genome.write_text('>chr1\nACGT\n')
            identity=root/'ref.json';identity.write_text(json.dumps(dict(reference_identity='ref',contigs=[
                dict(id='chr1',length=4,uppercase_sha256=hashlib.sha256(b'ACGT').hexdigest())])))
            bam=root/'test.bam';bam.write_bytes(b'fixture')
            request=encoded(dict(assembly_id='A',readset=dict(readset_id='R1',strandedness='forward',layout='single')))
            pysam=MagicMock();pysam.__version__='test';handle=pysam.AlignmentFile.return_value.__enter__.return_value
            handle.references=['chr1'];handle.lengths=[4];handle.header.to_dict.return_value={'HD':{'SO':'coordinate'}}
            with patch.dict(sys.modules,{'pysam':pysam}):
                check(genome,identity,bam,request,root/'checked.json')
                self.assertEqual(json.loads((root/'checked.json').read_text())['reference_identity'],'ref')
                self.assertEqual(json.loads((root/'checked.json').read_text())['orientation'],'SE')
                bad_orientation=encoded(dict(assembly_id='A',orientation='FR',readset=dict(readset_id='R1',strandedness='forward',layout='single')))
                with self.assertRaisesRegex(ValueError,'orientation does not match'):
                    check(genome,identity,bam,bad_orientation,root/'bad')
                handle.references=['other']
                with self.assertRaisesRegex(ValueError,'dictionary'):check(genome,identity,bam,request,root/'bad')
                genome.write_text('>chr1\nAAAA\n')
                with self.assertRaisesRegex(ValueError,'genome differs'):check(genome,identity,bam,request,root/'bad')

    def test_cluster_fixture_has_canonical_motifs_on_both_strands(self):
        from validate_junctions_cluster import fixtures
        from reference_repeats import fasta
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixtures(root)
            sequence=next(fasta(root/'accepted/repeats/A.masked.fasta'))[1]
            self.assertEqual(sequence[200:202]+sequence[298:300],'GTAG')
            self.assertEqual(sequence[1100:1102]+sequence[1198:1200],'CTAC')
            records=[s.split('\t') for s in (root/'sams/A.sam').read_text().splitlines() if not s.startswith('@')]
            self.assertEqual(len(records),252)
            self.assertTrue(all(len(r[9])==100 and len(r[10])==100 for r in records))

    def fixtures(self,root,strand='+',rna=True):
        # Coding junction 200..300 and UTR-associated junction 400..500 (half-open).
        gff=root/'A.gff3'
        records=[('gene',101,650,'.','ID=g1'),('mRNA',101,650,'.','ID=t1;Parent=g1'),
            ('exon',101,200,'.','ID=e1;Parent=t1'),('exon',301,400,'.','ID=e2;Parent=t1'),
            ('exon',501,650,'.','ID=e3;Parent=t1'),('CDS',151,200,'0','ID=c1;Parent=t1'),
            ('CDS',301,350,'1','ID=c2;Parent=t1')]
        gff.write_text('##gff-version 3\n'+''.join(f'chr1\tfixture\t{k}\t{s}\t{e}\t.\t{strand}\t{p}\t{a}\n' for k,s,e,p,a in records))
        identity=root/'identity.json';identity.write_text(json.dumps(dict(reference_identity='ref',contigs=[dict(id='chr1',length=1000)])))
        validation=root/'status.json';validation.write_text(json.dumps(dict(status='ready',stage='final_annotation',
            reference_identity='ref',gff3_sha256=checksum(gff))))
        readsets=[dict(readset_id=rid,layout='single',strandedness='reverse',library_id='L'+rid) for rid in ('R1','R2')] if rna else []
        request=encoded(dict(assembly_id='A',readsets=readsets))
        folders=[]
        for row in readsets:
            folder=root/row['readset_id'];folder.mkdir();folders.append(folder)
            (folder/'inputs.json').write_text(json.dumps(dict(status='ready',assembly_id='A',readset=row,reference_identity='ref')))
            (folder/'version.txt').write_text('portcullis 1.2.4\n')
            (folder/'filter.log').write_text('Will apply a lenient rule-based filter instead.\n')
            passing=row['readset_id']=='R1'
            self.native(folder/'all.junctions.tab',[(200,299,strand,12),(400,499,strand,2)])
            self.native(folder/'filtered.pass.junctions.tab',[(200,299,strand,12)] if passing else [])
            self.native(folder/'filtered.fail.junctions.tab',[(400,499,strand,2)]+([] if passing else [(200,299,strand,12)]))
            finish(folder,.5)
        return gff,validation,identity,request,folders

    def native(self,path,rows):
        path.write_text('refname\tstart\tend\tconsensus-strand\tnb_raw_aln\n'+''.join(f'chr1\t{s}\t{e}\t{strand}\t{n}\n' for s,e,strand,n in rows))

    def test_both_strands_exact_coordinates_library_union_and_coding_split(self):
        for strand in ('+','-'):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);args=self.fixtures(root,strand)
                result=report(*args,root/'out')
                self.assertEqual(result['counts']['pass'],1);self.assertEqual(result['counts']['fail'],1)
                self.assertTrue(result['readsets'][0]['native_low_junction_fallback'])
                with (root/'out/junction_support.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
                self.assertEqual((rows[0]['intron_start'],rows[0]['intron_end']),('201','300'))
                self.assertEqual(rows[0]['passing_readsets'],'R1');self.assertEqual(rows[0]['failed_readsets'],'R2')
                self.assertEqual(rows[0]['region'],'coding');self.assertEqual(rows[1]['region'],'utr_associated')
                self.assertEqual(rows[0]['raw_alignments'],'24')

    def test_unknown_or_opposite_strand_is_never_counted_as_support(self):
        for strand,label in (('.','unresolved_strand'),('-','not_observed')):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);g,v,i,r,folders=self.fixtures(root)
                for folder in folders:
                    self.native(folder/'all.junctions.tab',[(200,299,strand,12)])
                    self.native(folder/'filtered.pass.junctions.tab',[(200,299,strand,12)])
                    self.native(folder/'filtered.fail.junctions.tab',[]);finish(folder,.5)
                result=report(g,v,i,r,folders,root/'out')
                self.assertEqual(result['counts']['pass'],0)
                self.assertEqual(result['counts'][label],1 if strand=='.' else 2)

    def test_no_rna_is_not_assessed_and_never_zero_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);args=self.fixtures(root,rna=False);result=report(*args,root/'out')
            self.assertFalse(result['assessed']);self.assertEqual(result['counts']['not_assessed'],2)
            with (root/'out/transcript_junction_support.tsv').open() as h:row=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(row['all_pass_fraction'],'.');self.assertEqual(row['all_pass'],'.')

    def test_changed_tables_models_missing_readsets_and_partition_errors_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,v,i,r,folders=self.fixtures(root)
            with self.assertRaisesRegex(ValueError,'Missing/unexpected'):report(g,v,i,r,folders[:1],root/'bad')
            path=folders[0]/'filtered.pass.junctions.tab';original=path.read_text();path.write_text(original+'\n')
            with self.assertRaisesRegex(ValueError,'table changed'):report(g,v,i,r,folders,root/'bad')
            path.write_text(original)
            self.native(folders[0]/'filtered.fail.junctions.tab',[])
            with self.assertRaisesRegex(ValueError,'partition'):finish(folders[0],.5)
            g.write_text(g.read_text()+'# changed\n')
            with self.assertRaisesRegex(ValueError,'unchanged accepted'):report(g,v,i,r,folders,root/'bad')

    def test_one_base_shift_empty_tables_and_out_of_bounds_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,v,i,r,folders=self.fixtures(root)
            for folder in folders:
                self.native(folder/'all.junctions.tab',[(201,299,'+',12)])
                self.native(folder/'filtered.pass.junctions.tab',[(201,299,'+',12)])
                self.native(folder/'filtered.fail.junctions.tab',[]);finish(folder,.5)
            self.assertEqual(report(g,v,i,r,folders,root/'out')['counts']['not_observed'],2)
            path=root/'empty.tab';self.native(path,[]);self.assertEqual(table(path),{})
            self.native(path,[(900,1000,'+',1)])
            with self.assertRaisesRegex(ValueError,'outside accepted'):table(path,{'chr1':1000})

    def test_single_exon_has_no_fraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);g,v,i,r,folders=self.fixtures(root)
            rows=g.read_text().splitlines();g.write_text('\n'.join(x for x in rows if 'ID=e2;' not in x and 'ID=e3;' not in x and 'ID=c2;' not in x)+'\n')
            state=json.loads(v.read_text());state['gff3_sha256']=checksum(g);v.write_text(json.dumps(state))
            report(g,v,i,r,folders,root/'out')
            with (root/'out/transcript_junction_support.tsv').open() as h:row=next(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(row['status'],'single_exon');self.assertEqual(row['all_pass_fraction'],'.')


if __name__=='__main__':unittest.main()
