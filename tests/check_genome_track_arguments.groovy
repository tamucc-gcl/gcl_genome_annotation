// Run with GroovyMain and the Nextflow runtime on the Java classpath.
// Exercise the production command builders with actual Nextflow TaskPath values.
import nextflow.processor.TaskPath
import java.nio.file.Paths

def module = new File('modules/genomeTracks.nf').text
def commandSetup = { processName ->
    def process = module.substring(module.indexOf("process ${processName} {"))
    def start = process.indexOf('    script:') + '    script:'.length()
    process.substring(start, process.indexOf('"""', start))
}
def first = new TaskPath(Paths.get('native/E.out'), 'aux/E.out')
def second = new TaskPath(Paths.get('native/F.out'), 'aux/F.out')
// Reproduce the original failure: iterating a single path yields its components.
assert first.collect { it.toString() } == ['aux', 'E.out']
for (files in [[], first, [first], [first, second]]) {
    def expected = files instanceof Collection ? files : [files]
    def quoted = expected.collect { "'${it}'" }.join(' ')
    def filterBinding = new Binding([aux:files, policy:[enabled:true]])
    def extra = new GroovyShell(this.class.classLoader, filterBinding).evaluate(
        commandSetup('FILTER_GENOME_TRACK') + '\nextra')
    assert extra == quoted : "Filter arguments split staged paths: ${extra}"
    def reportBinding = new Binding([statuses:files, settings:[:]])
    def arguments = new GroovyShell(this.class.classLoader, reportBinding).evaluate(
        commandSetup('REPORT_GENOME_TRACKS') + '\narguments')
    assert arguments == quoted : "Report arguments split staged paths: ${arguments}"
}
println 'PASS: production filter/report arguments preserve single staged TaskPaths, empty lists and multiple files'
