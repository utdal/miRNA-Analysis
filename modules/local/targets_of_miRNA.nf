process TARGETS_OF_MIRNA {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'   

    input:
    path(miRNA_DE_file)

    output:
    tuple path("*all_up_reg_miRNA_targets_miRTarBase_TargetScan.tsv"), path("*all_down_reg_miRNA_targets_miRTarBase_TargetScan.tsv"),  emit: miRNA_targets
    tuple path("*up_regulated_miRNAs.tsv"), path("*down_regulated_miRNAs.tsv"),                                        optional: true, emit: deseq2_output_split
    path("*.log"),                                                                                                                     emit: target_database_log
    path "versions.yml",                                                                                                               emit: versions

    script:
    def hsa_miRTarBase_TargetScan_db = params.target_database
    def experimental_evidence = params.experimental_evidence ?: "weak"
    def min_weighted_context_percentile = params.min_weighted_context_percentile ?: 50
    def max_predicted_KD = params.max_predicted_KD ?: 0.0
    """
    get_miRNA_targets_from_database.py \
        --miRNA_list ${miRNA_DE_file} \
        --deseq2_output ${params.deseq2_output} \
        --hsa_miRTarBase_TargetScan_db ${hsa_miRTarBase_TargetScan_db} \
        --experimental_evidence ${experimental_evidence} \
        --min_weighted_context_percentile ${min_weighted_context_percentile} \
        --max_predicted_KD ${max_predicted_KD}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}