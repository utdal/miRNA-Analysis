//
// main workflow
//

// Import Local Subworkflows
include { FASTQC_CUTADAPT_FASTQC as PREPROCESSING} from './../subworkflows/local/preprocessing.nf'
//include { FASTQ_FIND_MIRNA_MIRDEEP2 } from './subworkflows/nf-core/fastq_find_mirna_mirdeep2/main'

// Import Local Modules
include { MIRBASE_DOWNLOAD } from './../modules/local/miRBase_download.nf'
include { EXCERPT } from './../modules/local/exceRpt.nf'
include { EXCERPTMERGE } from './../modules/local/exceRptMerge.nf'
include { HTSEQ_COUNT } from './../modules/local/htseq-count.nf'
include { CONCAT_RAW_COUNTS } from './../modules/local/concat_raw_counts.nf'
include { DESEQ2 } from './../modules/local/deseq2.nf'


workflow MIRNA_EXPRESSION {

    take:
    samplesheet
    meta_data

    main:

    ch_versions = Channel.empty()

    // Get input reads
    ch_samplesheet = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
    
    ch_samples = ch_samplesheet.map { row -> 
        def meta = [ id: row.sample_id, single_end: true ] 
        tuple(meta, file(row.fastq)) 
    }

    // preprocess: fastqc cutadapt fastqc
    PREPROCESSING (ch_samples)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())

    // Download miRBase files
    MIRBASE_DOWNLOAD ()
    ch_versions = ch_versions.mix(MIRBASE_DOWNLOAD.out.versions.first())

    //
    // miRDeep2
    //
    //FASTQ_FIND_MIRNA_MIRDEEP2 (
        //PREPROCESSING.out.trimmed_reads, 
        // genome fasta
        // bowtie index???
        // channel: [ val(meta),  mature_mirna, hairpin_mirna ]
    //)

    // exceRpt
    EXCERPT (
        PREPROCESSING.out.trimmed_reads
    )
    ch_versions = ch_versions.mix(EXCERPT.out.versions.first())

    // TODO Combine exceRpt results
    //EXCERPTMERGE (
        //EXCERPT.out.exceRpt_folder
    //)

    // HTSeq-count
    HTSEQ_COUNT (
        EXCERPT.out.exceRpt_aligned_bam,
        MIRBASE_DOWNLOAD.out.miRNA_gff
    )
    ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    
    // Concatenate raw counts
    CONCAT_RAW_COUNTS (
        HTSEQ_COUNT.out.raw_counts
    )
    ch_versions = ch_versions.mix(CONCAT_RAW_COUNTS.out.versions.first())

    // DESeq2???
    DESEQ2 (
        CONCAT_RAW_COUNTS.out.all_raw_counts,
        meta_data
    )

    emit:
    miRNA_DE = DESEQ2.out.miRNA_DE
    
}