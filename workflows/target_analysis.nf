
// Import Local Subworkflows
include { TARGETS } from './../subworkflows/local/targets.nf'

workflow TARGET_ANALYSIS {

    take:
    bulk_rna_counts
    miRNA_DE

    main:

    ch_versions = Channel.empty()

    TARGETS (
        miRNA_DE,
        bulk_rna_counts
    )
    ch_versions = ch_versions.mix(TARGETS.out.versions.first())


    // TODO miRanda or TargetScan for miRDeep2 results



    // Functional Analysis

    // PANTHER API: https://pantherdb.org/services/openAPISpec.jsp
    // PANTHER statistical enrichment using the Mann-Whitney U Test (Wilcoxon Rank-Sum Test).
    // For each term in the specified annotation data set, the genes associated with that term 
    //          are evaluated according to the likelihood that their numerical values were drawn 
    //          randomly from the overall distribution of values. 
    // In Linux sample curl command to invoke service is as follows: curl -X POST 'https://pantherdb.org/services/oai/pantherdb/enrich/statenrich' -H 'accept: application/json' -H 'Content-Type: multipart/form-data' -F 'organism=9606' -F 'annotDataSet=GO:0008150' -F 'correction=FDR' -F 'geneExp=@/path/to/two_column/tab/delimited/text_file.txt;type=text/plain' 
    // It is recommended that response from previous web service request is received before sending
    //           a new request. Failure to comply with this policy may result in the IP address 
    //           being blocked from accessing PANTHER.
}