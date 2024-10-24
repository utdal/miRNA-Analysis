process EXCERPT{
    tag "$meta.id"
    label 'process_medium'

    // The container 
    container "docker://rkitchen/excerpt"

    // TODO Get the reference files
    // make user download the files
    input:
    tuple val(meta), path(reads)
    path(ref_genome)

    // what was the input?????
    path()

    output:


    script:
    """
    
    
    
    """




}