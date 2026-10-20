// GroovyMain + Nextflow runtime: production channel routing and command construction.
import nextflow.Channel
import nextflow.processor.TaskPath
import nextflow.script.ScriptParser
import nextflow.script.ScriptBinding
import java.nio.file.Paths
import static nextflow.Nextflow.tuple
new nextflow.Session([:])
new ScriptParser(this.class.classLoader).setBinding(new ScriptBinding()).parse(new File('tests/release_fixture.nf').text)
def workflowText = new File('workflows/annotationRelease.nf').text
def routeCode = workflowText.substring(workflowText.indexOf('    optionalRows='),workflowText.indexOf('    releases='))
def refs=['B','A'].collect { a -> tuple(a,Paths.get(a+'.fa'),Paths.get(a+'.masked.fa'),Paths.get(a+'.json'),Paths.get(a+'.repeat.gff3'),Paths.get(a+'.repeat.json')) }
def preamble='import static nextflow.Nextflow.tuple\ndef file(Map opts,String path) { java.nio.file.Paths.get(path) }\n'
[true,false].each { enabled ->
    def captured=[]
    def binding=new Binding([
        params:[run_functional_annotation:enabled],projectDir:'.',
        references:Channel.fromList(refs),
        identified:Channel.fromList(['A','B'].collect { tuple(it,Paths.get(it+'/identified')) }),
        trackStatuses:Channel.fromList(['B','A'].collect { tuple(it,Paths.get(it+'.tracks.json')) }),
        tracks:Channel.fromList([tuple('A','telomere',Paths.get('A.telomere.gff3'))]),
        functional:Channel.fromList(enabled ? ['A','B'].collect { a -> tuple(a,*['coding','table','genes','status'].collect {Paths.get(a+'.'+it)}) } : []),
        PACKAGE_ANNOTATION_RELEASE:{ jobs, helpers -> captured.addAll(jobs.toList().getVal()) }
    ])
    new GroovyShell(this.class.classLoader,binding).evaluate(preamble+routeCode)
    assert captured.size()==2
    captured.each { r ->
        def a=r[0]
        assert r.size()==14 && r[1]==Paths.get(a+'/identified') && r[2].toString()==a+'.fa' && r[7].toString()==a+'.tracks.json' : r
        assert r[8].size()==(a=='A'?1:0) && r[9]==enabled : r
        assert enabled ? r[10].toString()==a+'.coding' : r[10]==[] : r
    }
}
// A single staged path must remain one auxiliary file, not its path components.
def module=new File('modules/annotationRelease.nf').text
def setup=module.substring(module.indexOf('    def trackFiles='),module.indexOf('    """'))
def path={ name -> new TaskPath(Paths.get(name),name) }
[true,false].each { enabled ->
    def vars=[aid:'A',identified:path('identified'),genome:path('genome.fa'),masked:path('masked.fa'),identity:path('identity.json'),
        repeats:path('repeat_source.gff3'),repeat_status:path('repeat.json'),track_status:path('tracks.json'),
        tracks:enabled ? path('selected_tracks/A.telomere.gff3') : [],functional_enabled:enabled,
        functional_gff:path('functional/coding.gff3'),functional_table:path('functional/table.tsv'),
        functional_genes:path('functional/genes.tsv'),functional_status:path('functional/status.json')]
    def result=new GroovyShell(this.class.classLoader,new Binding(vars)).evaluate(setup+'\nrequest')
    assert result.tracks==(enabled ? ['selected_tracks/A.telomere.gff3'] : [])
    assert result.functional==(enabled ? [gff:'functional/coding.gff3',table:'functional/table.tsv',genes:'functional/genes.tsv',status:'functional/status.json'] : [:])
}
println 'PASS: release fixture syntax, two-assembly native routing, functional on/off, empty/single-file staged paths'
System.exit(0)
