import base64
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from test_annotation_release import release_fixture
from package_annotation import package
from unified_report import build
from reference_repeats import sha


def report_fixture(root, mode='rna_protein', sources=None, functional=True):
    sources=sources or [release_fixture(root)]
    stages=[]; globals={k:dict(status='ready',assemblies=[]) for k in ('preflight','protein','rna','rna_alignment','reference_repeats','braker')}
    globals['rna']['readsets']=[];globals['rna_alignment']['readsets']=[]
    releases=[];structural=[]
    def save(key, aid, value):
        target=root/'statuses'/f'{key}.{aid}.json';target.parent.mkdir(exist_ok=True)
        target.write_text(json.dumps(value));stages.append([key,aid,str(target)]);return target
    for source in sources:
        aid=source['aid'];source=copy.deepcopy(source);ident=Path(source['identified']);identity=json.loads(Path(source['identity']).read_text());ref=identity['reference_identity']
        selected=ident.parent/'selection/selection_status.json';sel=json.loads(selected.read_text())
        if sel['assembly_id']!=aid:
            sel['assembly_id']=aid;selected.write_text(json.dumps(sel))
            a=ident/'annotation_status.json';state=json.loads(a.read_text());state['selection_sha256']=sha(selected);a.write_text(json.dumps(state))
        if source['functional_enabled']:
            p=Path(source['functional']['status']);f=json.loads(p.read_text());f.update(named_transcripts=1,hypothetical_transcripts=1,tools={
                'diamond':dict(hit=1,no_hit=1,ineligible=0,disabled=0), 'interproscan':dict(hit=0,no_hit=0,ineligible=0,disabled=2)})
            p.write_text(json.dumps(f))
        if not functional:source.update(functional_enabled=False,functional={})
        release=root/'release'/aid;package(source,release);releases.append(str(release))
        sq=root/'structural'/aid;sq.mkdir(parents=True);structural.append([aid,str(sq)])
        (sq/'completeness_metrics.tsv').write_text('model_set\tmetric\tvalue\nbefore\tbusco_complete\t90.0\nafter\tbusco_complete\t90.0\n')
        (sq/'plots').mkdir()
        png=base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+aZ1sAAAAASUVORK5CYII=')
        for name in ('retention','completeness','aed','junctions','structure','complexity','evidence_agreement','junction_occurrences','flags','decisions'):
            (sq/f'plots/{name}.png').write_bytes(png)
        s=dict(status='ready',stage='structural_qc',assembly_id=aid,reference_identity=ref,
            identified_model_sha256=sha(ident/f'{aid}.gff3'),model_sets=[dict(model_set='before',genes=2,transcripts=3),dict(model_set='after',genes=1,transcripts=2)],
            outputs={p.relative_to(sq).as_posix():sha(p) for p in sq.rglob('*') if p.is_file()})
        (sq/'structural_qc_status.json').write_text(json.dumps(s))
        assembly_mode=mode[aid] if isinstance(mode,dict) else mode
        has_rna=assembly_mode in ('rna','rna_protein');has_protein=assembly_mode in ('protein','rna_protein');readsets=['R1'] if has_rna else []
        globals['preflight']['assemblies'].append(dict(assembly_id=aid,settings=dict(species='Test species',taxid='1',telomere_motif='TTAGGG',busco_lineage='test_odb10')))
        globals['protein']['assemblies'].append(dict(assembly_id=aid,status='ready' if has_protein else 'not_requested',protein_count=12 if has_protein else 0))
        for key in ('rna','rna_alignment'):globals[key]['assemblies'].append(dict(assembly_id=aid,status='ready' if has_rna else 'not_requested',rna_readsets=readsets))
        if has_rna:
            globals['rna']['readsets'].append(dict(readset_id='R1',read_units=100))
            globals['rna_alignment']['readsets'].append(dict(assembly_id=aid,readset_id='R1',mapped_fraction=.75,strand=dict(assembly_id=aid,strandedness='reverse')))
        globals['reference_repeats']['assemblies'].append(dict(assembly_id=aid,status='ready',reference=dict(reference_identity=ref,stats=dict(excluded_contigs=0))))
        globals['braker']['assemblies'].append(dict(assembly_id=aid,status='ready',reference_identity=ref,mode=assembly_mode,genes=2,coding_transcripts=3))
        # Preserve the exact accepted selection manifest, not a reconstructed threshold table.
        stages.append(['selection',aid,str(ident.parent/'selection/selection_status.json')])
        save('utr',aid,dict(status='ready',reference_identity=ref,skipped=not has_rna,changed_transcripts=0,rejected_utr_candidates={},cds_sequences_unchanged=True,protein_sequences_unchanged=True))
    for key,state in globals.items():save(key,'run',state)
    return dict(releases=releases,structural=structural,stages=stages,configuration=dict(rfam_version='CURRENT',run_aed=True,run_junction_qc=False))


@unittest.skipUnless(importlib.util.find_spec('matplotlib'),'Unified report plotting environment required')
class UnifiedReportTests(unittest.TestCase):
    def test_four_modes_results_methods_only_and_assets(self):
        for mode in ('rna_protein','rna','protein','ab_initio'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);request=report_fixture(root,mode);out=root/'report';state=build(request,out)
                text=(out/'annotation_report.md').read_text(encoding='utf-8')
                self.assertIn('## Assembly 1: A',text);self.assertIn(mode,text)
                self.assertIn('| Metric | Before selection | Final selected |',text)
                self.assertIn('| BUSCO complete (%) | 90 | 90 |',text)
                self.assertIn('| AED | no | Not set | Not applied |',text)
                self.assertIn('| Coding junction support | no | Not set | Not applied |',text)
                self.assertIn('| tRNA | tRNA loci |',text)
                self.assertNotIn('| single_exon_short_cds |',text)
                self.assertNotIn('Tuning Recommendations',text);self.assertNotIn('consider tightening',text)
                self.assertIn('Rfam',text);self.assertIn('acquisition setting',text)
                self.assertEqual(state['assemblies'],['A'])
                self.assertFalse(list(out.rglob('*.html')))
                self.assertTrue((out/'annotation_report_plots/A.functional_coverage.svg').is_file())
                self.assertIn('not requested',text) if mode in ('protein','ab_initio') else self.assertIn('75.0%',text)

    def test_changed_release_and_missing_stage_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);r=report_fixture(root);bad=copy.deepcopy(r);bad['stages']=[v for v in bad['stages'] if v[0]!='rna']
            with self.assertRaisesRegex(ValueError,'Missing report stage'):build(bad,root/'bad')
            p=Path(r['releases'][0])/'tables/transcript_assessment.tsv';p.write_text(p.read_text()+'changed')
            with self.assertRaisesRegex(ValueError,'input changed'):build(r,root/'changed')

    def test_functional_disabled_and_no_rna(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);r=report_fixture(root,'ab_initio',functional=False);build(r,root/'report')
            text=(root/'report/annotation_report.md').read_text(encoding='utf-8')
            self.assertIn('Functional annotation: disabled.',text)
            self.assertFalse((root/'report/annotation_report_plots/A.functional_coverage.png').exists())


if __name__=='__main__':unittest.main()
