process CLUSTERPROFILER {
    label 'process_low'

    container "docker://utdpaincenter/mirna-analysis-cluterprofiler:1.0"

    input:
    tuple path(targets_up), path(targets_down)

    output:
    tuple path("all_up*"), path("all_down*"), emit: clusterprofiler_results
    path "clusterProfiler.log",               emit: clusterprofiler_log
    path "versions.yml",                      emit: versions

    script:
    def baseDir = params.base_dir
    """
    Rscript ${baseDir}/bin/clusterProfiler.R \\
        --up_miRNA_targets ${targets_up} \\
        --down_miRNA_targets ${targets_down}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$( R --version | awk 'NR==1{print \$3}' )
    END_VERSIONS
    """
}