process EXCERPTMERGE {
    label 'process_medium'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'

    input:
    val(all_folders)

    output:
    path("miRNA_counts.tsv"),               emit: miRNA_counts
    path("miRNA_rpm.tsv"),                  emit: miRNA_rpm
    path("piRNA_counts.tsv"),               emit: piRNA_counts
    path("piRNA_rpm.tsv"),                  emit: piRNA_rpm
    path("circRNA_counts.tsv"),             emit: circRNA_counts
    path("circRNA_rpm.tsv"),                emit: circRNA_rpm
    path("tRNA_counts.tsv"),                emit: tRNA_counts
    path("tRNA_rpm.tsv"),                   emit: tRNA_rpm
    path("GENCODE_counts.tsv"),             emit: GENCODE_counts
    path("GENCODE_rpm.tsv"),                emit: GENCODE_rpm
    path("rRNA_filtered_read_counts.tsv"),  emit: rRNA_filtered_read_counts
    path("biotype_counts.tsv"),             emit: biotype_counts
    path("merge_exceRpt_output.log"),       emit: excerpt_merged_log
    path "versions.yml",                    emit: versions

    script:
    """
    merge_exceRpt_output.py \\
        --exceRpt_results_dirs ${all_folders.join(',')}

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}