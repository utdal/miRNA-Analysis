//
// main workflow
//

// Import Local Subworkflows
include { FASTQC_CUTADAPT_FASTQC as PREPROCESSING} from './../subworkflows/local/preprocessing.nf'
include { FASTQ_FIND_MIRNA_MIRDEEP2 } from './../subworkflows/nf-core/fastq_find_mirna_mirdeep2/main.nf'

// Import Local Modules
include { MIRBASE_DOWNLOAD } from './../modules/local/miRBase_download.nf'
include { EXCERPT } from './../modules/local/exceRpt.nf'
include { EXCERPTMERGE } from './../modules/local/exceRptMerge.nf'
include { HTSEQ_COUNT } from './../modules/local/htseq-count.nf'
include { CONCAT_RAW_COUNTS } from './../modules/local/concat_raw_counts.nf'
include { DESEQ2 } from './../modules/local/deseq2.nf'
include { MERGEMIRDEEP2 } from '../modules/local/mergemiRDeep2.nf'


workflow MIRNA_EXPRESSION {

    take:
    skip_preprocessing
    skip_mirdeep2
    genome_fasta
    bowtie_index
    samplesheet
    meta_data_files

    main:

    ch_versions = Channel.empty()

    // TODO Add a path check for each input file
    // Get input reads
    ch_samplesheet = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
    
    ch_samples = ch_samplesheet.map { row -> 
        def meta = [ id: row.sample_id, single_end: true ] 
        tuple(meta, file(row.fastq)) 
    }

    if (skip_preprocessing != true) {
        // preprocess: fastqc cutadapt fastqc
        PREPROCESSING (ch_samples)
        reads = PREPROCESSING.out.trimmed_reads
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())
    } else {
        reads = ch_samples
    }
    

    // Download miRBase files
    MIRBASE_DOWNLOAD ()
    ch_versions = ch_versions.mix(MIRBASE_DOWNLOAD.out.versions.first())

    //
    // miRDeep2
    //
    if (skip_mirdeep2 != true) {
        // Make a channel for miRNA reference files
        ch_mature_hsa = MIRBASE_DOWNLOAD.out.mature_hsa_miRNA.map { it -> it}
        ch_hairpin_hsa = MIRBASE_DOWNLOAD.out.hairpin_hsa_miRNA.map { it -> it}
        ch_mature_other = MIRBASE_DOWNLOAD.out.mature_other_miRNA.map { it -> it}
        ch_mature_hairpin = ch_mature_hsa
            .combine(ch_hairpin_hsa)
            .combine(ch_mature_other)
            .map { mature_hsa, hairpin_hsa, mature_other -> 
                [[id: 'mature_hairpin'], mature_hsa, hairpin_hsa, mature_other] 
            }
            .first()
        
        // Run miRDeep2 and then merge the results
        FASTQ_FIND_MIRNA_MIRDEEP2 (
            reads, 
            genome_fasta,
            bowtie_index,
            ch_mature_hairpin
        )
        ch_versions = ch_versions.mix(FASTQ_FIND_MIRNA_MIRDEEP2.out.versions)

        // Get all miRDeep2 results output files
        all_miRDeep2_results = FASTQ_FIND_MIRNA_MIRDEEP2.out.outputs.collect{it[1]}
        
        // Merge miRDeep2 results
        MERGEMIRDEEP2 (
            all_miRDeep2_results
        )
        ch_versions = ch_versions.mix(MERGEMIRDEEP2.out.versions.first())
    }

    // exceRpt
    EXCERPT (
        reads
    )
    ch_versions = ch_versions.mix(EXCERPT.out.versions.first())

    all_excerpt_folders = EXCERPT.out.exceRpt_folder.collect{it[1]}
    EXCERPTMERGE (
        all_excerpt_folders
    )
    

    // HTSeq-count then merge results
    HTSEQ_COUNT (
        EXCERPT.out.exceRpt_aligned_bam,
        MIRBASE_DOWNLOAD.out.miRNA_gff
    )
    | combine
    | collect
    | CONCAT_RAW_COUNTS
    ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    ch_versions = ch_versions.mix(CONCAT_RAW_COUNTS.out.versions.first())

    // Parse all of the meta datas
    ch_meta_data_files = Channel
        .fromPath(meta_data_files)
        .splitCsv(header: true)
    
    ch_meta_data = ch_meta_data_files.map { row -> 
        def meta2 = [ condition: row.condition ] 
        tuple(meta2, file(row.csv)) 
    }

    // DESeq2
    DESEQ2 (
        CONCAT_RAW_COUNTS.out.all_raw_counts,
        ch_meta_data
    )

    emit:
    miRNA_DE = DESEQ2.out.miRNA_DE  // channel: [ val(meta2), path("*.tsv") ]
    versions = ch_versions
    
}