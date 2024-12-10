include { MIRNA_EXPRESSION } from './workflows/miRNA_Expression'
include { TARGET_ANALYSIS } from './workflows/target_analysis'

workflow MIRNA_ANALYSIS {

    // TODO make optional
    MIRNA_EXPRESSION(
        params.samplesheet,
        params.meta_data_files
    )  

    // TODO make optional
    TARGET_ANALYSIS(
        MIRNA_EXPRESSION.out.miRNA_DE,
        params.bulk_rna_counts
    )
}
