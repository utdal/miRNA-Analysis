include { MIRNA_EXPRESSION } from './workflows/miRNA_Expression'
include { TARGET_ANALYSIS } from './workflows/target_analysis'

workflow {

    versions = Channel.empty()

    // Prep mirdeep2 inputs
    genome_fasta = params.genome_fasta? Channel
                                            .fromPath(params.genome_fasta, checkIfExists: true)
                                            .map{ it -> [ [id:it.baseName], it ] }
                                            .collect() : Channel.empty()
    bowtie_index = params.bowtie_index ? Channel
                                            .value( [ [id:""], file(params.bowtie_index) ] ) : Channel.empty()
    
    // TODO make optional
    MIRNA_EXPRESSION(
        params.skip_preprocessing,
        params.skip_mirdeep2,
        genome_fasta,
        bowtie_index,
        params.samplesheet,
        params.meta_data_files
    )
    versions = versions.mix(MIRNA_EXPRESSION.out.versions)

    // TODO make optional
    TARGET_ANALYSIS(
        params.miRNA_DE,
        params.bulk_rna_counts
    )
    versions = versions.mix(TARGET_ANALYSIS.out.versions)
}
