process DESEQ2 {
    label 'process_medium'

    // Container obtained from htseq-count module from nf-core
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'biocontainers/bioconductor-deseq2:1.34.0--r41hc247a5b_3' }"

    input:
    path(all_raw_counts)
    path(metadata)

    output:
    path("*.tsv"), emit: intersect_targets

    script:
    """
    #!/bin/Rscript
    library(DESeq2)
    # TODO add the rest of the commands

    """


}