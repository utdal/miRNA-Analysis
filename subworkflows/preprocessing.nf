//
// QC Trim QC
//

include { FASTQC as FASTQC_RAW  } from './../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from './../modules/nf-core/fastqc/main'
include { CUTADAPT as CUTADAPT_ADAPTER } from './../modules/nf-core/cutadapt/main'
include { CUTADAPT as CUTADAPT_EXTRA } from './../modules/nf-core/cutadapt/main'

workflow FASTQC_CUTADAPT_FASTQC{

    ch_versions = Channel.empty()

    take:
    ch_reads            // channel: [ val(meta), path(reads)]

    main:

    // Raw quality
    // TODO make optional
    FASTQC_RAW ( ch_reads )
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

    emit:
    trimmed_reads    = CUTADAPT_EXTRA.out.reads // channel: [ val(meta), path(html) ]
    
    fastqc_raw_html  = FASTQC_RAW.out.html      // channel: [ val(meta), path(html) ]
    fastqc_raw_zip   = FASTQC_RAW.out.zip       // channel: [ val(meta), path(zip) ]
    fastqc_trim_html = FASTQC_TRIM.out.html     // channel: [ val(meta), path(html) ]
    fastqc_trim_zip  = FASTQC_TRIM.out.zip      // channel: [ val(meta), path(zip) ]

    versions   = ch_versions.ifEmpty(null)      // channel: [ path(versions.yml) ]
    
}