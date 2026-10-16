import hashlib
import base64
import json
import re
import shutil
import subprocess
import textwrap
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from utr_refinement import finalize, validate, dump, prepare, boundary_label_completions
from prepare_protein_evidence import checksum
from gene_model_aed import gene_scores
from evidence_aed import models
from inspect_utr_labels import inspect as inspect_labels


class UtrContracts(unittest.TestCase):
    def test_label_diagnostic_identifies_cds_at_terminal_exon_base(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original,raw,refined,identity,genome,inputs=self.fixture(root)
            def write(path,rows):
                path.write_text(''.join(f'chr1\tx\t{k}\t{s}\t{e}\t.\t+\t{phase}\t{attrs}\n' for k,s,e,phase,attrs in rows))
            coding=[('CDS',101,101,'0','ID=c1;Parent=t1'),('CDS',201,299,'2','ID=c2;Parent=t1')]
            write(original,[('gene',101,299,'.','ID=g1'),('mRNA',101,299,'.','ID=t1;Parent=g1'),
                ('exon',101,101,'.','ID=e1;Parent=t1'),('exon',201,299,'.','ID=e2;Parent=t1')]+coding)
            write(raw,[('gene',81,320,'.','ID=g1'),('mRNA',81,320,'.','ID=t1;Parent=g1;utr_refine_evidence=r1'),
                ('exon',81,101,'.','ID=e1;Parent=t1'),('exon',201,320,'.','ID=e2;Parent=t1'),
                ('three_prime_UTR',300,320,'.','ID=u3;Parent=t1')]+coding)
            result=inspect_labels(original,raw,identity,root/'labels.json')
            self.assertEqual(result['mismatched_labels_on_eligible_candidates'],1)
            self.assertEqual(result['issues'][0]['missing_label_blocks'],[(80,100)])
            self.assertTrue(result['issues'][0]['boundary_equalities']['cds_start_is_exon_last_base'])
            state=json.loads(inputs.read_text());state['model_sha256']=checksum(original);dump(inputs,state)
            finalize(original,raw,identity,refined)
            status=self.run_validation(root,(original,raw,refined,identity,genome,inputs))
            self.assertEqual(status['completed_boundary_utr_labels'],{'t1':{'five_prime_UTR':[[80,100]]}})
            self.assertEqual(models(refined,{'chr1':1000})['t1']['exons'],[(80,101),(200,320)])

    def test_boundary_label_completion_is_strand_correct_and_narrow(self):
        for strand in ('+','-'):
            low_kind='five_prime_UTR' if strand=='+' else 'three_prime_UTR'
            high_kind='three_prime_UTR' if strand=='+' else 'five_prime_UTR'
            observed={low_kind:[],high_kind:[(250,270)]}
            result=boundary_label_completions([(80,101),(200,270)],[(100,101),(200,250)],observed,[],strand)
            self.assertEqual(result,{low_kind:[(80,100)]})
            observed={low_kind:[(80,100)],high_kind:[]}
            result=boundary_label_completions([(80,150),(200,230)],[(100,150),(200,201)],observed,[(201,204)],strand)
            self.assertEqual(result,{high_kind:[(204,230)]})
            # A missing ordinary UTR label is not this known boundary bug.
            self.assertEqual(boundary_label_completions([(80,150),(200,270)],[(100,150),(200,250)],
                {low_kind:[],high_kind:[]},[],strand),{})
            # Extra unexpected UTR labels remain a validation failure.
            self.assertEqual(boundary_label_completions([(80,101),(200,270)],[(100,101),(200,250)],
                {low_kind:[(70,75)],high_kind:[(250,270)]},[],strand),{})

    def test_label_diagnostic_excludes_separate_stop_codons(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original,raw,refined,identity,genome,inputs=self.fixture(root)
            self.assertEqual(inspect_labels(original,raw,identity,root/'labels.json')['mismatched_labels'],0)

    def test_native_bam_list_command_preserves_multiple_rows_without_heredoc(self):
        module=(Path(__file__).resolve().parents[1]/'modules/utrRefinement.nf').read_text()
        body=module.split('process RUN_UTR_REFINEMENT {',1)[1].split('process FINALIZE_UTR_ANNOTATION',1)[0]
        command=re.search(r'python -c "(import base64;[^\n]+)"',body).group(1)
        for count in (0,1,2,10):
            with self.subTest(readsets=count),tempfile.TemporaryDirectory() as tmp:
                expected=''.join(f'bams/readset{i:02}.bam\ttrue\ttrue\n' for i in range(1,count+1))
                encoded=base64.b64encode(expected.encode()).decode()
                result=subprocess.run([sys.executable,'-c',command.replace('${fofEncoded}',encoded)],cwd=tmp,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertEqual((Path(tmp)/'bam.fof').read_bytes(),expected.encode())

    def test_preparation_rejects_unresolved_bam_protocol_and_reports_rna_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); original,raw,refined,identity,genome,inputs=self.fixture(root)
            validation=root/'validation.json'; dump(validation,dict(status='ready',gff3_sha256=checksum(original),reference_identity='ref'))
            metadata=root/'readsets.json';dump(metadata,[dict(layout='paired',strandedness='reverse')])
            bam=root/'reads.bam'; bam.touch(); index=root/'reads.bam.csi';index.touch()
            rna=root/'rna.gtf';rna.write_text(''.join(f'chr1\tStringTie\t{k}\t81\t220\t.\t{strand}\t.\tgene_id "{tid}"; transcript_id "{tid}";\n' for tid,strand in [('r1','+'),('r2','.')] for k in ('transcript','exon')))
            prepare(original,validation,identity,rna,metadata,[bam],[index],root/'prep')
            state=json.loads((root/'prep/inputs.json').read_text())
            self.assertEqual(state['rna_transcripts'],1)
            self.assertEqual(state['rna_ambiguous_strand_excluded'],1)
            dump(metadata,[dict(layout='paired',strandedness='forward')])
            with self.assertRaisesRegex(ValueError,'adapter'): prepare(original,validation,identity,rna,metadata,[bam],[index],root/'prep')

    def fixture(self,root,strand='+'):
        genome=root/'genome.fa'; genome.write_text('>chr1\n'+'A'*1000+'\n')
        identity=root/'identity.json'; dump(identity,dict(reference_identity='ref',contigs=[dict(id='chr1',length=1000,uppercase_sha256=hashlib.sha256(b'A'*1000).hexdigest())]))
        original=root/'original.gff3'
        rows=[('gene',101,202,'.','ID=g1'),('mRNA',101,202,'.','ID=t1;Parent=g1'),
              ('exon',101,199,'.','ID=e1;Parent=t1'),('CDS',101,199,'0','ID=c1;Parent=t1')]
        if strand=='+': rows.append(('stop_codon',200,202,'0','ID=stop1;Parent=t1'))
        def write(path,rows):
            path.write_text('##gff-version 3\n'+''.join(f'chr1\tx\t{k}\t{s}\t{e}\t.\t{strand}\t{p}\t{a}\n' for k,s,e,p,a in rows))
        write(original,rows)
        raw=root/'native.gff3'
        write(raw,[('gene',81,220,'.','ID=g1'),('mRNA',81,220,'.','ID=t1;Parent=g1;utr_refine_evidence=r1'),
            ('exon',81,220,'.','ID=ne1;Parent=t1'),('CDS',101,199,'0','ID=c1;Parent=t1'),
            ('five_prime_UTR' if strand=='+' else 'three_prime_UTR',81,100,'.','ID=u1;Parent=t1'),
            ('three_prime_UTR' if strand=='+' else 'five_prime_UTR',200,220,'.','ID=u2;Parent=t1')])
        inputs=root/'inputs.json'; dump(inputs,dict(model_sha256=checksum(original),reference_identity='ref',readsets=[]))
        for name in ('old.cds','new.cds','old.proteins','new.proteins'): (root/name).write_text('>t1\nAAA\n')
        refined=root/'refined.gff3'; finalize(original,raw,identity,refined)
        return original,raw,refined,identity,genome,inputs

    def run_validation(self,root,files,skipped=False):
        original,raw,refined,identity,genome,inputs=files
        return validate(original,refined,identity,genome,inputs,root/'old.cds',root/'new.cds',root/'old.proteins',root/'new.proteins',root/'out',skipped,refined.with_suffix('.rejections.json'))

    def test_wrong_strand_native_proposal_restores_only_affected_isoform(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); files=self.fixture(root)
            original,raw,refined,identity,genome,inputs=files
            def second_isoform(text,wrong_strand=False):
                result=[]
                for line in text.splitlines():
                    f=line.split('\t')
                    if len(f)!=9 or f[2]=='gene': continue
                    f[8]=f[8].replace('t1','t2').replace('ne1','ne2').replace('ID=e1','ID=e2').replace('ID=c1','ID=c2').replace('ID=u1','ID=u3').replace('ID=u2','ID=u4').replace('stop1','stop2')
                    if wrong_strand and f[2]=='exon': f[6]='-'
                    result.append('\t'.join(f))
                return '\n'.join(result)+'\n'
            original.write_text(original.read_text()+second_isoform(original.read_text()))
            raw.write_text(raw.read_text()+second_isoform(raw.read_text(),True))
            state=json.loads(inputs.read_text());state['model_sha256']=checksum(original);dump(inputs,state)
            for name in ('old.cds','new.cds','old.proteins','new.proteins'):
                with (root/name).open('a') as h:h.write('>t2\nAAA\n')
            finalize(original,raw,identity,refined)
            before=models(original,{'chr1':1000});after=models(refined,{'chr1':1000})
            self.assertEqual(after['t2']['exons'],before['t2']['exons'])
            self.assertNotEqual(after['t1']['exons'],before['t1']['exons'])
            status=self.run_validation(root,files)
            self.assertEqual(set(status['rejected_utr_candidates']),{'t2'})
            self.assertEqual(status['changed_transcripts'],1)
            self.assertEqual(after['t2']['strand'],'+')

    def test_native_extensions_preserve_cds_and_exclude_separate_stop(self):
        for strand in ('+','-'):
            with self.subTest(strand=strand),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); files=self.fixture(root,strand); status=self.run_validation(root,files)
                self.assertEqual(status['changed_transcripts'],1)
                import csv
                with (root/'out/utr_changes.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(row['five_prime_utr_bp'],'20' if strand=='+' else '21')
                self.assertEqual(row['three_prime_utr_bp'],'18' if strand=='+' else '20')
                if strand=='+': self.assertIn('stop_codon',files[2].read_text())

    def test_unmarked_native_change_restores_original_and_reports_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);files=self.fixture(root)
            original,raw,refined,identity,genome,inputs=files
            raw.write_text(raw.read_text().replace(';utr_refine_evidence=r1',''))
            finalize(original,raw,identity,refined)
            status=self.run_validation(root,files)
            self.assertEqual(status['changed_transcripts'],0)
            self.assertIn('lacks UTR evidence provenance',status['rejected_utr_candidates']['t1'])
            self.assertEqual(models(original,{'chr1':1000})['t1']['exons'],models(refined,{'chr1':1000})['t1']['exons'])

    def test_unmarked_duplicate_exon_rows_restore_original_without_counting_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);files=self.fixture(root)
            original,raw,refined,identity,genome,inputs=files
            duplicate=next(line for line in original.read_text().splitlines() if '\texon\t' in line).replace('ID=e1','ID=duplicate')
            raw.write_text(original.read_text()+duplicate+'\n')
            finalize(original,raw,identity,refined)
            status=self.run_validation(root,files)
            self.assertEqual(status['changed_transcripts'],0)
            self.assertEqual(set(status['rejected_utr_candidates']),{'t1'})

    def test_restoring_unmarked_candidate_cannot_hide_native_cds_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original,raw,refined,identity,genome,inputs=self.fixture(root)
            raw.write_text(raw.read_text().replace(';utr_refine_evidence=r1','').replace('\t+\t0\tID=c1','\t+\t1\tID=c1'))
            with self.assertRaisesRegex(ValueError,'Native UTR refinement changed CDS'):
                finalize(original,raw,identity,refined)

    def test_no_rna_passthrough_and_sequence_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); files=self.fixture(root); files[2].write_bytes(files[0].read_bytes())
            status=self.run_validation(root,files,True)
            self.assertTrue(status['skipped']); self.assertEqual(status['changed_transcripts'],0)
            (root/'new.proteins').write_text('>t1\nAAC\n')
            with self.assertRaisesRegex(ValueError,'sequences'): self.run_validation(root,files,True)

    def test_cds_phase_and_reference_changes_block_refinement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); files=self.fixture(root)
            files[2].write_text(files[2].read_text().replace('\t+\t0\tID=c1','\t+\t1\tID=c1'))
            with self.assertRaisesRegex(ValueError,'phases'): self.run_validation(root,files)
            finalize(files[0],files[1],files[3],files[2])
            files[4].write_text('>chr1\n'+'C'*1000+'\n')
            with self.assertRaisesRegex(ValueError,'genome'): self.run_validation(root,files)


class OptionalAedFilterContracts(unittest.TestCase):
    def test_existing_filter_keeps_default_models_and_uses_selected_source(self):
        source=Path(__file__).resolve().parents[1]
        template=(source/'modules/filterGeneModels.nf').read_text().split('    """')[1]
        for selected,expected in [(None,{'g1','g2'}),('rna',{'g2'}),('protein',{'g1'})]:
            with self.subTest(source=selected),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)
                gff=root/'models.gff3'
                gff.write_text(''.join(f'chr1\tx\t{k}\t{s}\t{e}\t.\t+\t.\t{attrs}\n'
                    for gid,tid,s,e in [('g1','t1',101,500),('g2','t2',701,1100)]
                    for k,attrs in [('gene','ID='+gid),('mRNA','ID='+tid+';Parent='+gid),
                        ('exon','ID='+tid+'.e;Parent='+tid),('CDS','ID='+tid+'.c;Parent='+tid)]))
                proteins=root/'proteins.fa';proteins.write_text('>t1\nMAAA\n>t2\nMCCC\n')
                scores=root/'scores.tsv';scores.write_text('gene_id\ttranscript_id\trna_aed\tprotein_aed\ng1\tt1\t0.9\t0\ng2\tt2\t0.1\t0.9\n')
                values=dict(assembly_id='A',min_cds_length='0',min_gene_length='0',min_intron_length='0',
                    braker_gff3=gff.as_posix(),braker_proteins=proteins.as_posix(),aed_scores=scores.as_posix())
                values["params.filter_max_aed != null ? params.filter_max_aed : 'None'"]='0.7' if selected else 'None'
                values['params.filter_aed_source ?: ""']=selected or ''
                script=textwrap.dedent(re.sub(r'\$\{([^}]+)\}',lambda m:values[m[1]],template)).replace('\\\\','\\')
                (root/'filter.py').write_text(script)
                shutil.copy2(source/'py_scripts/gene_model_aed.py',root/'gene_model_aed.py')
                result=subprocess.run([sys.executable,str(root/'filter.py')],cwd=root,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                import csv
                with (root/'A.filtered_aed.tsv').open() as handle:
                    reader=csv.DictReader(handle,delimiter='\t');rows=list(reader)
                    self.assertEqual(reader.fieldnames,['gene_id','transcript_id','rna_aed','protein_aed'])
                self.assertEqual({r['gene_id'] for r in rows},expected)

    def test_disabled_filter_does_not_require_evidence(self):
        self.assertEqual(gene_scores('does-not-exist',None,None,{'g1'}),{})

    def test_sources_are_separate_and_missing_source_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scores.tsv'; p.write_text('gene_id\ttranscript_id\trna_aed\tprotein_aed\ng1\tt1\t0.9\t0\ng1\tt2\t0.5\t0.2\n')
            self.assertEqual(gene_scores(p,'rna',.7,{'g1'}),{'g1':.5})
            self.assertEqual(gene_scores(p,'protein',0,{'g1'}),{'g1':0})
            with self.assertRaisesRegex(ValueError,'filter_aed_source'): gene_scores(p,None,.7,{'g1'})
            p.write_text(p.read_text().replace('\t0.9\t','\t.\t'))
            with self.assertRaisesRegex(ValueError,'unavailable'): gene_scores(p,'rna',.7,{'g1'})


if __name__=='__main__': unittest.main()
