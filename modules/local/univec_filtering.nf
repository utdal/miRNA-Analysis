process UNIVEC_FILTERING {
    tag "$meta.id"
    label 'process_high'



    input:
    tuple val(meta), path(trimmed_reads)

    output:
    tuple val(meta), path("${meta.id}_univec_filtered/${meta.id}_univec_filtered.fastq.gz"), emit: univec_filtered_reads
    tuple val(meta), path("${meta.id}_univec_filtered"),          emit: univec_filtered_folder


    script:
    base_dir = params.base_dir
    align_threads = params.align_threads
    """

    make -f ${base_dir}/bin/exceRpt_univec_filtering \
        OUTPUT_DIR=${meta.id}_univec_filtered \
        SAMPLE_ID=${meta.id} \
        INPUT_FILE_PATH=${trimmed_reads} \
        BASE_DIR=${base_dir} \
        N_THREADS=${align_threads}
    
    
    """

}