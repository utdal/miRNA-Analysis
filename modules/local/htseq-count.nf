process HTSEQ_COUNT{
    tag "$meta.id"
    label 'process_medium'

    // Container obtained from htseq-count module from nf-core
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.3--py310ha14a713_0':
        'biocontainers/htseq:2.0.3--py310ha14a713_0' }"

    input:
    tuple val(meta), path(bam)
    path(gff3)

    output:
    tuple val(meta), path("*.tsv"), emit: raw_counts
    path "versions.yml",            emit: versions

    script:
    def args = task.ext.args ?: ''
    def threads = params.htseq_threads
    """
    htseq-count \\
        -f bam \\
        -t miRNA \\
        -i Name \\
        ${args} \\
        -c ${meta.id}_raw_counts.tsv \\
        -n ${threads} \\
        ${bam} \\
        ${gff3}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq-count: \$( htseq-count --version | sed 's/htseq-count\s\+//g' )        
    END_VERSIONS
    """

}