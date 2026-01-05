process UNIVEC_FILTERING {
    tag "$meta.id"
    label 'process_high'

    container "docker://utdpaincenter/mirna-analysis-excerpt-univec-filtering:1.2"

    input:
    tuple val(meta), path(trimmed_reads)

    output:
    tuple val(meta), path("${meta.id}_univec_filtered/${meta.id}_univec_filtered.fastq.gz"), emit: univec_filtered_reads
    tuple val(meta), path("${meta.id}_univec_filtered"),                                     emit: univec_filtered_folder
    path "versions.yml",                                                                     emit: versions


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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        star: \$(STAR --version | sed -n 's/STAR_//p' )
        samtools: \$( samtools --version | awk 'NR==1{print \$2}' )
        bowtie2: \$( bowtie2 --version | grep 'version' | awk '{print \$3}' )
        java-jdk: \$( java -version 2>&1 | awk -F '"' '/version/ {print \$2}' )
        r: \$( R --version | awk 'NR==1{print \$3}' )
    END_VERSIONS
    """

}