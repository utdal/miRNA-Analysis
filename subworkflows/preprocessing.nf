//
// QC Trim QC
//

include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { CUTADAPT as CUTADAPT_ADAPTER } from '../../../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_EXTRA } from '../../../modules/nf-core/cutadapt/main'

workflow FASTQC_CUTADAPT_FASTQC{

    ch_versions = Channel.empty()

    take:
    ch_reads            // channel: [ val(meta), path(reads)]

    // Raw quality
    // TODO make optional
    FASTQC_RAW ( ch_reads )
    FASTQC_RAW.out.versions.view()
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    // Remove adapters
    CUTADAPT_ADAPTER ( ch_reads )
    ch_versions = ch_versions.mix(CUTADAPT_ADAPTER.out.versions.first())

    // Remove extra stuff
    // TODO make optional
    CUTADAPT_EXTRA ( CUTADAPT_ADAPTER.out.reads )
    ch_versions = ch_versions.mix(CUTADAPT_EXTRA.out.versions.first())

    // Trimmed quality
    FASTQC_TRIM ( CUTADAPT_EXTRA.out.reads )
    ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())

    
}