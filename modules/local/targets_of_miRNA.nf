process TARGETS_OF_MIRNA {
    tag "$meta2.condition"
    label 'process_low'    

    input:
    tuple val(meta2), path(miRNA_DE_file)

    output:
    tuple val(meta2), path("*all_up_regulated_ENCORI_miRNA_targets.tsv"), path("*all_down_regulated_ENCORI_miRNA_targets.tsv"),  emit: miRNA_targets
    tuple val(meta2), path("*up_regulated_miRNAs.tsv"), path("*down_regulated_miRNAs.tsv"), optional: true, emit: deseq2_output_split
    tuple val(meta2), path("*.log"),  emit: target_database_log
    path "versions.yml",              emit: versions

    script:
    // TODO change output file name so that it is based on the meta2.condition
    def args = task.ext.args ?: ''
    """    
    get_miRNA_targets_from_database.py \
        --miRNA_list ${miRNA_DE_file} \
        --condition ${meta2.condition} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}