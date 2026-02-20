/*
========================================================================================
    DOWNLOAD InterProScan DATABASE MODULE
========================================================================================
    Downloads the specified InterProScan lineage database once, before parallel InterProScan runs.
    
    The downloaded database is stored in the specified output directory and can be
    reused across pipeline runs via Nextflow's caching.
========================================================================================
*/

process DOWNLOAD_INTERPROSCAN_DB {
    tag "${interproscan_version}"
    label 'interproscan_download'
    
    storeDir "${params.interproscan_downloads}"
    
    input:
    val interproscan_version
    
    output:
    path("interproscan-${interproscan_version}/data"), emit: db
    
    script:
    """
    curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${interproscan_version}/alt/interproscan-data-${interproscan_version}.tar.gz
    tar -pxzf interproscan-data-${interproscan_version}.tar.gz
    
    echo "Successfully downloaded InterProScan"
    """
    
    stub:
    """
    mkdir -p interproscan-data-${interproscan_version}/data
    """
}