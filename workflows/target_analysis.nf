//
// miRNA Target Analysis Workflow
//

// Import Local Modules
include { TARGETS_OF_MIRNA } from './../modules/local/targets_of_miRNA'
include { INTERSECT_MIRNA_RNASEQ } from './../modules/local/intersect_miRNA_RNAseq'
include { CLUSTERPROFILER } from './../modules/local/clusterProfiler'

workflow TARGET_ANALYSIS {

    take:
    miRNA_DE
    bulk_rna_counts

    main:
    ch_versions = Channel.empty()

    // Get all the targets of the differentially expressed miRNAs
    TARGETS_OF_MIRNA( 
        miRNA_DE 
    )
    ch_versions = ch_versions.mix(TARGETS_OF_MIRNA.out.versions.first())

    // Filter all the targets for only those that appear in the Bulk RNA-seq 
    //      (aka only genes that are expressed in tissue)
    // If no bulk_rna_counts is provided, then just return the unique targets
    INTERSECT_MIRNA_RNASEQ( 
        TARGETS_OF_MIRNA.out.miRNA_targets,
        bulk_rna_counts
    )
    ch_versions = ch_versions.mix(INTERSECT_MIRNA_RNASEQ.out.versions.first())

    CLUSTERPROFILER(
        INTERSECT_MIRNA_RNASEQ.out.filtered_intersect_targets
    )
    ch_versions = ch_versions.mix(CLUSTERPROFILER.out.versions.first())

    // TODO miRanda or TargetScan for miRDeep2 results
    
    emit:
    versions = ch_versions

}