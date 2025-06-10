process DESEQ2 {
    tag "$meta2.condition"
    label 'process_low'

    container "docker://utdpaincenter/mirna-analysis-deseq2:1.1"

    input:
    path(all_raw_counts)
    tuple val(meta2), path(metadata)

    output:
    tuple val(meta2), path("*_deseq2_results/*"),  emit: miRNA_DE
    tuple val(meta2), path("*_basic_plots_*.pdf"), emit: basic_plots
    tuple val(meta2), path("*_R_sessionInfo.log"), emit: R_sessionInfo
    tuple val(meta2), path("*_DESeq2.log"),        emit: DESeq2_log
    path "versions.yml",                           emit: versions

    script:
    def baseDir = params.base_dir
    def min_padj = params.min_padj ?: 0.05
    def min_lfc = params.min_lfc ?: 0.58
    """
    Rscript ${baseDir}/bin/deseq2.R \\
        --all_raw_counts ${all_raw_counts} \\
        --metadata ${metadata} \\
        --min_padj ${min_padj} \\
        --min_log2fc ${min_lfc} \\
        --meta2_condition ${meta2.condition}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$( R --version | awk 'NR==1{print \$3}' )
        deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}