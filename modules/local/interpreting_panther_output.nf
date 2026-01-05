process INTERPRETING_PANTHER_OUTPUT {
    label "process_medium"

    container "docker://utdpaincenter/mirna-analysis-interpreting-panther:1.0"

    input:
    tuple path(up_json), path(down_json)

    output:
    path("*/*"),            emit: panther_results
    path "versions.yml",    emit: versions

    script:
    """
    interpreting_PANTHER_output.py \\
        --up_panther_output "${up_json}" \\
        --down_panther_output "${down_json}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}