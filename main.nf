//
// main workflow
//

// Import Subworkflows
include { FASTQC_CUTADAPT_FASTQC as PREPROCESSING} from './subworkflows/preprocessing'

// Import Modules



workflow {

    ch_versions = Channel.empty()

    // Get inputs

    config_directory = params.config_directory
    reads = params.reads


    // TODO Check inputs


    // preprocess fastqc cutadapt fastqc
    // TODO reads needs to be a channel. 
    PREPROCESSING (reads)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())

    // miRDeep2



    // exceRpt



    // HTSeq-count




    // DESeq2???



    // TODO miRanda or TargetScan for miRDeep2 results



    // Functional Analysis
    
}