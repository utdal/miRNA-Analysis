#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MIRNA_EXPRESSION } from './workflows/miRNA_Expression'
include { TARGET_ANALYSIS } from './workflows/target_analysis'

workflow {

    versions = Channel.empty()

    // Only running miRNA expression workflow
    if (params.workflow.equals("mirna_expression")) {

        skip_preprocessing = params.skip_preprocessing ?: false
        skip_mirdeep2 = params.skip_mirdeep2 ?: false
        skip_deseq2 = params.skip_deseq2 ?: false

        // Prep mirdeep2 inputs
        if (skip_mirdeep2 != true) {
            genome_fasta = params.genome_fasta? Channel
                                                    .fromPath(params.genome_fasta, checkIfExists: true)
                                                    .map{ it -> [ [id:it.baseName], it ] }
                                                    .collect() : Channel.empty()
            bowtie_index = params.bowtie_index ? Channel
                                                    .value( [ [id:""], file(params.bowtie_index) ] ) : Channel.empty()
        }

        MIRNA_EXPRESSION(
            skip_preprocessing,
            skip_mirdeep2,
            skip_deseq2,
            genome_fasta,
            bowtie_index,
            params.samplesheet,
            params.meta_data_files
        )
        versions = versions.mix(MIRNA_EXPRESSION.out.versions)
    }

    // Only running target analysis workflow
    if (params.workflow.equals("target_analysis")) {
        TARGET_ANALYSIS(
            params.miRNA_DE,
            params.bulk_rna_counts
        )
        versions = versions.mix(TARGET_ANALYSIS.out.versions)
    }
}
