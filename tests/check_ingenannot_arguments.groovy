// Native AED preparation: the five-field production tuple and staged path arguments.
import nextflow.processor.TaskPath
import nextflow.script.ScriptParser
import nextflow.script.ScriptBinding
import java.nio.file.Paths
new nextflow.Session([:])
def text=new File('modules/ingenannotAed.nf').text
new ScriptParser(this.class.classLoader).setBinding(new ScriptBinding()).parse(text)
assert text.contains("tuple val(aid), path(gff), path(validation), path(rna, stageAs: 'rna/*'), path(protein, stageAs: 'protein/*')")
def setup=text.substring(text.indexOf('    def r ='),text.indexOf('    """'))
def path={ name -> new TaskPath(Paths.get(name),name) }
[[[],[]],[path('rna/merged.gtf'),[]],[[],path('protein/alignments.gff')],
 [[path('rna/merged.gtf')],[path('protein/alignments.gff')]]].each { pair ->
    def result=new GroovyShell(this.class.classLoader,new Binding([rna:pair[0],protein:pair[1]])).evaluate(setup+'\n[r,p,options]')
    assert result[0].size()==(pair[0] ? 1 : 0)
    assert result[1].size()==(pair[1] ? 1 : 0)
    assert !result[2].contains('current-status')
}
assert text.contains('structural_annotation/aed/${aid}/native')
println 'PASS: native AED five-field tuple, no-RNA/protein modes and single/multiple staged-path arguments'
System.exit(0)
