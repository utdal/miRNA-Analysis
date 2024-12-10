process TARGETS_OF_MIRNA {
    tag "$meta2.condition"
    label 'process_medium'    

    input:
    tuple val(meta2), path(miRNA_DE)

    output:
    tuple val(meta2), path("*.tsv"),          emit: all_targets
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