//
//
//

include { GET_MIRNA_TARGETS } from './../../modules/local/get_miRNA_targets'
include { INTERSECT_MIRNA_RNASEQ } from './../../modules/local/intersect_miRNA_RNAseq'

workflow TARGETS{

    ch_versions = Channel.empty()

    take:
    miRNA_DE

    main:

    // Get all the targets of the differentially expressed miRNAs
    GET_MIRNA_TARGETS( miRNA_DE )
    ch_versions = ch_versions.mix(GET_MIRNA_TARGETS.out.versions.first())

    // Filter all the targets for only those that appear in the Bulk RNA-seq 
    //      (aka only the ones that are expressed)
    INTERSECT_MIRNA_RNASEQ( GET_MIRNA_TARGETS.out.all_targets)
    ch_versions = ch_versions.mix(INTERSECT_MIRNA_RNASEQ.out.versions.first())

    // TODO condense the intersect_targets so that we have a list of genes with the miRNAs they target
    // TODO condense intersect_targets so that we have a list of the miRNAs and all the genes they target
    //      aiming for a GMT format maybe????



    emit:
    all_miRNA_DE_targets = GET_MIRNA_TARGETS.out.all_targets
    intersected_miRNA_DE_targets = INTERSECT_MIRNA_RNASEQ.out.intersect_targets

    versions = ch_versions.ifEmpty(null)

}