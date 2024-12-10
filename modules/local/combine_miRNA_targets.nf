process COMBINE_MIRNA_TARGETS {
    label 'process_medium'

    input:
    path(miRNA_DE_targets)

    output:
    path("*.tsv"),          emit: all_targets

    script:
    """
    """
}