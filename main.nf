include { MIRNA_EXPRESSION } from './workflows/miRNA_Expression'
include { TARGET_ANALYSIS } from './workflows/target_analysis'

workflow MIRNA_EXPRESSION_RUN {
    MIRNA_EXPRESSION()  
}

workflow TARGET_ANALYSIS_RUN {
    TARGET_ANALYSIS()
}

workflow {
    MIRNA_EXPRESSION_RUN()
}