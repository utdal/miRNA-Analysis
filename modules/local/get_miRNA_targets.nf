process GET_MIRNA_TARGETS {
    label 'process_medium'

    // TODO how to handle up vs. down regulated stuff

    input:
    path(miRNA_DE)

    output:
    path(all_targets), emit: all_targets

    script:
    """    
    python3 get_miRNA_targets_from_database.py \
        --miRNA_list ${miRNA_DE} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """


}