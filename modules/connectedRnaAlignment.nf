// The integrated entry consumes native reference/RNA outputs, not published paths.
process PLAN_CONNECTED_RNA_ALIGNMENT {
    cache false
    input:
    path preflight
    path passed
    path prepared
    path rna_passed
    path references, stageAs: 'reference_inputs/*'
    val locations
    path helper
    output:
    path 'alignment_plan.json', emit: plan
    script:
    """
    python ${helper} plan --preflight '${preflight}/resolved_inputs.json' --prepared '${prepared}' --request '${locations}' --output alignment_plan.json
    """
}
