process TARGETS_OF_MIRNA {
    label 'process_medium'    

    input:
    path(miRNA_DE)

    output:
    path("*.tsv"),          emit: all_targets
    path "versions.yml",    emit: versions

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