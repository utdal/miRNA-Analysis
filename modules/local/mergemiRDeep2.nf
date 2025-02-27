process MERGEMIRDEEP2 {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'

    input:
    path(outputs)

    output:
    path "similar_hairpin_all_samples.csv",                      emit: similar_hairpin_all_samples
    path "only_one_novel_miRNA_detected_across_all_samples.csv", emit: only_one_novel_miRNA_detected_across_all_samples
    path "all_samples.csv",                                      emit: all_samples
    path "*.log",                                                emit: mergemirdeep2_log
    path "versions.yml",                                         emit: versions

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