process EXCERPT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_exceRpt"),                                                        emit: exceRpt_folder
    tuple val(meta), path("${meta.id}_exceRpt/${meta.id}/endogenousAlignments_genome_Aligned.out.bam"), emit: exceRpt_aligned_bam
    path "versions.yml",                                                                                emit: versions

    script:
    base_dir = params.base_dir
    align_threads = params.align_threads
    """
    make -f ${base_dir}/bin/exceRpt_smallRNA \
        OUTPUT_DIR=${meta.id}_exceRpt \
        SAMPLE_ID=${meta.id} \
        INPUT_FILE_PATH=${reads} \
        BASE_DIR=${base_dir} \
        N_THREADS=${align_threads} \
        JAVA_RAM=16g
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        star: \$(STAR --version | sed -n 's/STAR_//p' )
        samtools: \$( samtools --version | awk 'NR==1{print \$2}' )
        bowtie2: \$( bowtie2 --version | grep 'version' | awk '{print \$3}' )
        java-jdk: \$( java -version 2>&1 | awk -F '"' '/version/ {print \$2}' )
        r: \$( R --version | awk 'NR==1{print \$3}' )
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}