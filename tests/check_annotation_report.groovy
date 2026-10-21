// Run with GroovyMain and the Nextflow runtime: production report path arguments.
import nextflow.processor.TaskPath
import nextflow.script.ScriptParser
import nextflow.script.ScriptBinding
import java.nio.file.Paths

new ScriptParser(this.class.classLoader).setBinding(new ScriptBinding()).parse(new File('tests/annotation_report_fixture.nf').text)
def module=new File('modules/annotation_report.nf').text
def setup=module.substring(module.indexOf('    def asList'),module.indexOf('    """'))
def path={ p -> new TaskPath(Paths.get(p),p) }
[false,true].each { multiple ->
    def r=multiple ? [path('release_inputs/A'),path('release_inputs/B')] : path('release_inputs/A')
    def s=multiple ? [path('structural_inputs/A'),path('structural_inputs/B')] : path('structural_inputs/A')
    def stages=multiple ? [path('stage01/same.json'),path('stage02/same.json')] : path('stage01/same.json')
    def vars=[releases:r,structural:s,structural_ids:multiple?['A','B']:['A'],stages:stages,
              stage_keys:multiple?[['utr','A'],['utr','B']]:[['utr','A']],configuration:[:]]
    def request=new GroovyShell(this.class.classLoader,new Binding(vars)).evaluate(setup+'\nrequest')
    assert request.releases==(multiple?['release_inputs/A','release_inputs/B']:['release_inputs/A'])
    assert request.structural[0]==['A','structural_inputs/A']
    assert request.stages[0]==['utr','A','stage01/same.json']
    if (multiple) assert request.stages[1]==['utr','B','stage02/same.json']
}
println 'PASS: report fixture syntax and production single/multiple staged-path arguments, including repeated source filenames'
