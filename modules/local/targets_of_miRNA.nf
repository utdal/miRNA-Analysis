process TARGETS_OF_MIRNA {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'   

    input:
    path(miRNA_DE_file)

    output:
    tuple path("*all_up_regulated_ENCORI_miRNA_targets.tsv"), path("*all_down_regulated_ENCORI_miRNA_targets.tsv"),  emit: miRNA_targets
    tuple path("*up_regulated_miRNAs.tsv"), path("*down_regulated_miRNAs.tsv"), optional: true, emit: deseq2_output_split
    path("*.log"),  emit: target_database_log
    path "versions.yml",              emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    get_miRNA_targets_from_database.py \
        --miRNA_list ${miRNA_DE_file} \
        --deseq2_output ${params.deseq2_output} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}