process MIRBASE_DOWNLOAD {
    label 'process_low'

    output:
    path("hsa.gff3"),               emit: miRNA_gff
    path("mature_hsa.fa"),          emit: mature_hsa_miRNA
    path("hairpin_hsa.fa"),         emit: hairpin_hsa_miRNA
    path("mature_other.fa"),        emit: mature_other_miRNA
    path("get_miRBase_files.log"),  emit: miRBase_logger
    path "versions.yml",            emit: versions

    script:
    base_dir = params.base_dir
    """
    python ${base_dir}/bin/get_miRBase_files.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}