//
// main workflow
//

// Import Subworkflows
include { FASTQC_CUTADAPT_FASTQC as PREPROCESSING} from './subworkflows/preprocessing.nf'

// Import Modules



workflow {

    ch_versions = Channel.empty()

    // Get inputs
    //config_directory = params.config_directory
    
    // TODO Take input as a csv file or something???
    def meta = [
        id: 'trial',
        single_end: true
    ]
    values = tuple(meta, file(params.reads))

    // TODO Check inputs


    // preprocess fastqc cutadapt fastqc
    // TODO reads needs to be a channel. 
    PREPROCESSING (values)
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())

    // miRDeep2



    // exceRpt



    // HTSeq-count




    // DESeq2???



    // TODO miRanda or TargetScan for miRDeep2 results



    // Functional Analysis
    
}