// Run from the checkout with GroovyMain and the Nextflow runtime on the classpath.
import nextflow.script.ScriptParser
import nextflow.script.ScriptBinding

new nextflow.Session([:])
def root = new File('.').absoluteFile.toPath().normalize().toFile()
def entries = new File(root, 'validation/entrypoints').listFiles().findAll { it.name.endsWith('.nf') }
assert entries.size() == 13
assert root.listFiles().findAll { it.name.endsWith('.nf') }*.name == ['main.nf']

[root, new File(root, 'validation/entrypoints')].each { directory ->
    def parser = new ConfigSlurper()
    parser.setBinding([projectDir: directory.toPath(), launchDir: root.toPath()])
    def config = parser.parse(new File(root, 'nextflow.config').getText('UTF-8'))
    assert config.params.pipeline_root == root.toString()
    ['ingenannot', 'utr_ranking', 'protein_evidence_alignment', 'agat'].each { label ->
        def environment = config.process[label].conda.toString()
        assert environment.startsWith(root.toString()) : environment
        assert new File(environment).isFile() : environment
    }
}

def scripts = [new File(root, 'main.nf')] + entries +
    ['modules', 'workflows', 'tests'].collectMany { directory ->
        new File(root, directory).listFiles().findAll { it.name.endsWith('.nf') }
    }
scripts.each { script ->
    def text = script.getText('UTF-8')
    new ScriptParser(this.class.classLoader).setBinding(new ScriptBinding()).parse(text)
    (text =~ /from\s+['"]([^'"]+)['"]/).each { match ->
        // Fixture launchers are copied to a staged pipeline root by their controllers.
        def includeRoot = script.parentFile.name == 'tests' ? root : script.parentFile
        def target = new File(includeRoot, match[1])
        if (!target.name.endsWith('.nf')) target = new File(target.toString() + '.nf')
        assert target.isFile() : "${script}: ${target}"
    }
    // Literal staged helper paths must exist after resolving the shared checkout root.
    (text =~ /\$\{params\.pipeline_root\}\/([^'"\s$]+\.(?:py|json|pl))/).each { match ->
        assert new File(root, match[1]).isFile() : "${script}: ${match[1]}"
    }
}
println "PASS: ${scripts.size()} Nextflow scripts/modules, nested/root config paths, includes and helpers"
