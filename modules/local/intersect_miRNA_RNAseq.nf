process INTERSECT_MIRNA_RNASEQ{
    label 'process_low'

    input:
    path(all_targets)

    output:
    path(intersect_targets), emit: intersect_targets

    script:
    """
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}