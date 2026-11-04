// Run with the Nextflow runtime, like the other check_*.groovy checks.
import nextflow.script.ScriptParser
import nextflow.script.ScriptBinding
new nextflow.Session([:])

def aed = new File('modules/evidenceAed.nf').text.replace('\r\n', '\n')
def utr = new File('modules/utrRefinement.nf').text.replace('\r\n', '\n')
def junction = new File('workflows/junctionSupport.nf').text.replace('\r\n', '\n')
[aed, utr, new File('modules/rnaAlignment.nf').text].each {
    new ScriptParser(this.class.classLoader).setBinding(new ScriptBinding()).parse(it)
}
def assembly = aed.split(java.util.regex.Pattern.quote('process ASSEMBLE_AED_RNA {'), 2)[1].split('\n}', 2)[0].split('    script:\n', 2)[1]
def routing = junction.substring(junction.indexOf('            def strand ='), junction.indexOf('            def request ='))
def ranking = utr.split(java.util.regex.Pattern.quote('process RUN_UTR_REFINEMENT {'), 2)[1].split('    script:\n', 2)[1]
ranking = ranking.substring(0, ranking.indexOf('    def nativeCommand'))

['forward', 'reverse', 'unstranded'].each { setting ->
    def command = new GroovyShell(this.class.classLoader, new Binding([
        strand: setting, bam: 'r.bam', rid: 'R', task: [cpus: 2]
    ])).evaluate(assembly)
    assert command.contains('--fr') == (setting == 'forward')
    assert command.contains('--rf') == (setting == 'reverse')
}
// A fallback carries unresolved inference but a valid operational setting.
def row = [readset_id:'R', layout:'paired', inferred_strandedness:'unresolved',
           effective_strandedness:'unstranded', strandedness:'unstranded',
           effective_setting_source:'conservative_fallback']
['FR', 'RF'].each { geometry ->
    def values = new GroovyShell(this.class.classLoader, new Binding([
        row:row, params:[portcullis_paired_orientation:geometry]
    ])).evaluate(routing + '\n[strand,orientation]')
    assert values == ['unstranded', geometry]
}
def fof = new GroovyShell(this.class.classLoader, new Binding([
    readsets:[row], bams:['r.bam'], rna:[]
])).evaluate(ranking + '\nfofEncoded')
assert new String(fof.decodeBase64(), 'UTF-8') == 'r.bam\ttrue\tfalse\n'
println 'PASS: StringTie flags, unstranded Portcullis with independent mate geometry, InGenAnnot fallback FOF'
System.exit(0)
