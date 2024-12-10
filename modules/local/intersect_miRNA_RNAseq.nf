process INTERSECT_MIRNA_RNASEQ {
    tag "$meta2.condition"
    label 'process_low'

    input:
    tuple val(meta2), path(all_targets)
    path(tissue_specific_genes)

    output:
    tuple val(meta2), path("*.tsv"),       emit: intersect_targets
    path "versions.yml",                   emit: versions

    script:
    """
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}