/*
========================================================================================
    PARSE SAMPLE SHEET FUNCTION (CACHE-OPTIMIZED)
========================================================================================
    Reads and validates the input sample sheet CSV file
    
    Expected CSV columns:
        - assembly_id: Unique sample identifier
        - assembly_path: Path to genome assembly FASTA file
    
    Rows with missing information are skipped and logged
    
    CACHING BEHAVIOR:
    This function creates independent channel items per sample using .fromList(),
    which ensures that adding new samples to the samplesheet will NOT invalidate
    the cache for existing samples when using -resume.
    
    PATH RESOLUTION:
    Relative paths in the samplesheet are resolved against the samplesheet's
    parent directory. Absolute paths (starting with "/") are used as-is.
========================================================================================
*/

def parseSampleSheet(sample_sheet_path) {
    
    def samples = []
    def sample_sheet_file = file(sample_sheet_path)
    
    if (!sample_sheet_file.exists()) {
        error "Sample sheet does not exist: ${sample_sheet_path}"
    }
    
    // Resolve relative paths against the samplesheet's parent directory
    def sheet_dir = sample_sheet_file.parent
    log.info "Resolving samplesheet paths relative to: ${sheet_dir}"

    def resolvePath = { String p ->
        if (!p) return null
        return p.startsWith('/') ? file(p) : file("${sheet_dir}/${p}")
    }
    
    def lines = sample_sheet_file.readLines()
    
    lines.eachWithIndex { line, idx ->
        // Skip header
        if (idx == 0 || line.startsWith('assembly_id')) {
            log.info "Parsing sample sheet: ${sample_sheet_path}"
            return
        }
        
        // Skip empty lines
        if (!line.trim()) return
        
        // Parse CSV line
        def fields = line.split(',').collect { it.trim() }
        
        if (fields.size() != 2) {
            log.warn "Line ${idx + 1}: Malformed CSV line (expected 2 fields, got ${fields.size()}): ${line}"
            return
        }
        
        def (assembly_id, assembly_path) = fields
        
        // Check for missing required fields
        def missing = []
        if (!assembly_id) missing << "assembly_id"
        if (!assembly_path) missing << "assembly_path"
        
        if (missing.size() > 0) {
            log.warn "Line ${idx + 1}: Skipping sample '${assembly_id ?: 'UNKNOWN'}': missing required field(s): ${missing.join(', ')}"
            return
        }
        
        // Resolve relative to samplesheet directory
        def assembly_file = resolvePath(assembly_path)
        
        if (!assembly_file.exists()) {
            log.warn "Line ${idx + 1}: Skipping sample '${assembly_id}': assembly file does not exist: ${assembly_path} (resolved to: ${assembly_file})"
            return
        }
        
        samples << tuple(assembly_id, assembly_file)
        log.info "Line ${idx + 1}: Successfully parsed assembly '${assembly_id}'"
    }
    
    if (samples.size() == 0) {
        error "No valid samples found in sample sheet: ${sample_sheet_path}"
    }
    
    log.info "Successfully parsed ${samples.size()} assembly(ies) from sample sheet"
    
    // CRITICAL: Use Channel.fromList() for independent cache entries per sample
    return Channel.fromList(samples)
}