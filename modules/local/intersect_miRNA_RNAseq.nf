process INTERSECT_MIRNA_RNASEQ {
    tag "$meta2.condition"
    label 'process_low'

    container 'docker://utdpaincenter/mirna-analysis-pandas-biopython-requests:1.3'
    

    input:
    tuple val(meta2), path(up_regulated_miRNA_targets), path(down_regulated_miRNA_targets)
    path(tissue_specific_genes)

    output:
    tuple val(meta2), path("*_up_regulated_final_target_gene_targeted_count.tsv"),  path("*_down_regulated_final_target_gene_targeted_count.tsv"), emit: intersect_targets
    tuple val(meta2), path("*.log"),       emit: intersect_targets_log
    path "versions.yml",                   emit: versions

    script:
    def condition = meta2.condition
    def args = task.ext.args ?: ''
    """

    intersect_miRNA_targets_with_BulkRNAseq.py \
        --up_miRNA_targets ${up_regulated_miRNA_targets} \
        --down_miRNA_targets ${down_regulated_miRNA_targets} \
        --bulkRNAseq ${tissue_specific_genes} \
        --condition ${condition} \
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}