process MERGEMIRDEEP2 {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'

    input:
    path(outputs)

    output:
    path "similarSort_hairpin_novel_miRNA.csv",           emit: similar_sort_hairpin
    path "similarSort_only_one_novel_miRNA_detected.csv", emit: similar_sort_only_one_novel_miRNA_detected_across
    path "exactSort_hairpin_mature_novel_miRNA.csv",      emit: exact_sort_novel_miRNAs
    path "exactSort_only_one_novel_miRNA_detected.csv",   emit: exact_sort_only_one_novel_miRNA_detected
    path "all_samples_novel_miRNAs.csv",                  emit: all_samples
    path "*.log",                                         emit: mergemirdeep2_log
    path "versions.yml",                                  emit: versions

    script:
    def min_mirdeep2_score = params.min_mirdeep2_score ?: 4.0
    def list_of_results = outputs.join(',')

    """
    merge_mirdeep2_results.py \\
        --miRDeep2_results_list ${list_of_results} \\
        --min_miRDeep2_score ${min_mirdeep2_score}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}