/*




*/

process EXCERPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${modulesDir}/local/exceRpt.yml"
    
    // Could not get it to work with singularity
    //container = 'docker://rkitchen/excerpt'

    // TODO Get the reference files
    // make user download the files
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_exceRpt")

    script:
    base_dir = params.base_dir

    align_threads = params.align_threads
    """

    make -f ${base_dir}/bin/exceRpt_smallRNA \
        OUTPUT_DIR=./${meta.id}_exceRpt \
        SAMPLE_ID=${meta.id} \
        READS=${reads} \
        BASE_DIR=${base_dir} \
        N_THREADS=${align_threads} \
        JAVA_RAM=16g

    """
    }