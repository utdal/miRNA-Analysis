process DESEQ2 {
    tag "$meta2.condition"
    label 'process_medium'

    // Container obtained from htseq-count module from nf-core
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'biocontainers/bioconductor-deseq2:1.34.0--r41hc247a5b_3' }"

    input:
    path(all_raw_counts)
    tuple val(meta2), path(metadata)

    output:
    tuple val(meta2), path("*.tsv"), emit: miRNA_DE

    script:
    """
    # TODO Call the script with DESeq2 commands
    #!/usr/bin/env R
    library(DESeq2)



    """


}