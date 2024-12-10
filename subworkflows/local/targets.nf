//
// miRNA Targets
//

// Import Local Modules
include { TARGETS_OF_MIRNA } from './../../modules/local/targets_of_miRNA'
include { INTERSECT_MIRNA_RNASEQ } from './../../modules/local/intersect_miRNA_RNAseq'

workflow TARGETS{
    take:
    miRNA_DE
    tissue_specific_genes

    main:
    ch_versions = Channel.empty()
    // Get all the targets of the differentially expressed miRNAs
    TARGETS_OF_MIRNA( miRNA_DE )
    ch_versions = ch_versions.mix(TARGETS_OF_MIRNA.out.versions.first())

    // Filter all the targets for only those that appear in the Bulk RNA-seq 
    //      (aka only the ones that are expressed)
    INTERSECT_MIRNA_RNASEQ( 
        TARGETS_OF_MIRNA.out.all_targets,
        tissue_specific_genes
    )
    ch_versions = ch_versions.mix(INTERSECT_MIRNA_RNASEQ.out.versions.first())

    // TODO condense the intersect_targets so that we have a list of genes with the miRNAs they target
    // TODO include a count of the up and down regulated miRNAs that target the gene
    // TODO condense intersect_targets so that we have a list of the miRNAs and all the genes they target
    //      aiming for a GMT format maybe????
    
    emit:
    all_miRNA_DE_targets = TARGETS_OF_MIRNA.out.all_targets
    intersected_miRNA_DE_targets = INTERSECT_MIRNA_RNASEQ.out.intersect_targets

    versions = ch_versions.ifEmpty(null)

}