process MIRBASE_DOWNLOAD {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'

    output:
    path("hsa.gff3"),               emit: miRNA_gff
    path("mature_hsa.fa"),          emit: mature_hsa_miRNA
    path("hairpin_hsa.fa"),         emit: hairpin_hsa_miRNA
    path("mature_other.fa"),        emit: mature_other_miRNA
    path("get_miRBase_files.log"),  emit: miRBase_logger
    path "versions.yml",            emit: versions

    script:
    """
   get_miRBase_files.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}