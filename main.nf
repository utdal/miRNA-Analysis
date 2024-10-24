//
// main workflow
//

// Import Subworkflows
include { FASTQC_CUTADAPT_FASTQC as PREPROCESSING} from './subworkflows/local/preprocessing.nf'
//include { FASTQ_FIND_MIRNA_MIRDEEP2 }              from './subworkflows/nf-core/fastq_find_mirna_mirdeep2/main'

// Import Modules


workflow {

    ch_versions = Channel.empty()

    // TODO Take input as a csv file or something???
    def meta = [
        id: 'trial',
        single_end: true
    ]
    reads = tuple(meta, file(params.reads))

    // TODO Check inputs


    // preprocess: fastqc cutadapt fastqc
    // TODO reads needs to be a channel. 
    PREPROCESSING (reads)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())

    
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



    // HTSeq-count




    // DESeq2???



    // TODO miRanda or TargetScan for miRDeep2 results



    // Functional Analysis
    
}