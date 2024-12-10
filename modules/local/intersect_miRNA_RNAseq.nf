process INTERSECT_MIRNA_RNASEQ{
    label 'process_low'

    input:
    path(all_targets)
    path(tissue_specific_genes)

    output:
    path("*.tsv"),          emit: intersect_targets
    path "versions.yml",    emit: versions

    script:
    """
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}