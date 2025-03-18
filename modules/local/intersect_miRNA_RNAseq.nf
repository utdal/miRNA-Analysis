process INTERSECT_MIRNA_RNASEQ {
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'
    

    input:
    tuple path(up_regulated_miRNA_targets), path(down_regulated_miRNA_targets)
    val tissue_specific_genes

    output:
    tuple path("up_regulated_final_target_gene_targeted_count.tsv"),  path("down_regulated_final_target_gene_targeted_count.tsv"),                                  emit: intersect_targets
    tuple path("RNAseq_filtered_up_miRNA_DE_targets.tsv"), path("RNAseq_filtered_down_miRNA_DE_targets.tsv"), optional: true,                                       emit: RNAseq_filtered_targets
    tuple path("filtered_up_regulated_final_target_gene_targeted_count.tsv"), path("filtered_down_regulated_final_target_gene_targeted_count.tsv"), optional: true, emit: filtered_intersect_targets
    path("*.log"),                                                                                                                                                  emit: intersect_targets_log
    path "versions.yml",                                                                                                                                            emit: versions

    script:
    def min_expression = params.min_expression ?: 0
    """

    intersect_miRNA_targets_with_BulkRNAseq.py \
        --up_miRNA_targets ${up_regulated_miRNA_targets} \
        --down_miRNA_targets ${down_regulated_miRNA_targets} \
        --bulkRNAseq ${tissue_specific_genes} \
        --min_expression ${min_expression}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}