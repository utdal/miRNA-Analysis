process PANTHER {
    label "process_low"

    input:
    tuple val(meta2), file(targets_up), file(targets_down)

    output:
    tuple val(meta2), path("*_up_results.json"), path("*_down_results.json"), emit: panther_results

    script:
    def condition = meta2.condition

    // Organism will be hardcoded to be for humans
    def organism = "9606"

    // annotDataSet options:
    // GO:0003674 - Molecular Function
    // GO:0008150 - Biological Process
    // GO:0005575 - Cellular Component
    // ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF - PANTHER GO Slim Molecular Function
    // ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP - PANTHER GO Slim Biological Process
    // ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC - PANTHER GO Slim Cellular Component
    // ANNOT_TYPE_ID_PANTHER_PC - protein class
    // ANNOT_TYPE_ID_PANTHER_PATHWAY - ANNOT_TYPE_PANTHER_PATHWAY
    // ANNOT_TYPE_ID_REACTOME_PATHWAY - ANNOT_TYPE_REACTOME_PATHWAY
    def annotDataSet = params.panther_annotDataSet ?: "GO:0008150"

    // correction options:
    // FDR - False Discovery Rate
    // BONFERRONI - Bonferroni
    // None - No correction
    def correction = params.panther_correction ?: "FDR"

    // A maximum of 100,000 genes can be submitted at a time
    """
    #!/bin/bash

    # Do input checking
    valid_annotations=("GO:0008150" "GO:0005575" "GO:0003674" \\
              "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF" "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" \\
              "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC" "ANNOT_TYPE_ID_PANTHER_PC" \\
              "ANNOT_TYPE_ID_PANTHER_PATHWAY" "ANNOT_TYPE_ID_REACTOME_PATHWAY")
              
    if [[ " \${valid_annotations[@]} " =~ " ${annotDataSet} " ]]; then
        echo "Valid annotation data set"
    else
        echo "Invalid annotation data set. Please choose from the following:"
        for value in "\${valid_annotations[@]}"; do
            echo "\$value"
        done
        exit 1
    fi

    if [[ "${correction}" == "FDR" || "${correction}" == "BONFERRONI" || "${correction}" == "None" ]]; then
        echo "Valid correction set."
    else
        echo "Invalid correction"
        exit 1
    fi

    # Check that the targets_up have less than 100,000 genes
    if [ \$(wc -l < "${targets_up}") -gt 100000 ]; then
        echo "Up regulated miRNA target genes file has more than 100,000 genes"
        echo "Please take the file and run directly in PANTHER website."
        skip_up=1
    fi

    # Check that the targets_down have less than 100,000 genes
    if [ \$(wc -l < ${targets_down}) -gt 100000 ]; then
        echo "Down regulated miRNA target genes file has more than 100,000 genes"
        echo "Please take the file and run directly in PANTHER website."
        skip_down=1
    fi

    echo "Input checks have been passed."

    # Check if the up reglated file is empty
    if [[ -s "${targets_up}" && "\$skip_down" -eq 0 ]]; then
        curl -X POST 'https://pantherdb.org/services/oai/pantherdb/enrich/statenrich' \\
            -H 'accept: application/json' \\
            -H 'Content-Type: multipart/form-data' \\
            -F 'organism=${organism}' \\
            -F 'annotDataSet=${annotDataSet}' \\
            -F 'correction=${correction}' \\
            -F 'geneExp=@${targets_up};type=text/plain' > ${condition}_up_results.json
    else
        echo "No up regulated miRNA target genes to perform analysis on."
    fi


    # Check if the down regulated file is empty
    if [[ -s "${targets_down}" && "\$skip_down" -eq 0 ]]; then
        curl -X POST 'https://pantherdb.org/services/oai/pantherdb/enrich/statenrich' \\
            -H 'accept: application/json' \\
            -H 'Content-Type: multipart/form-data' \\
            -F 'organism=${organism}' \\
            -F 'annotDataSet=${annotDataSet}' \\
            -F 'correction=${correction}' \\
            -F 'geneExp=@${targets_down};type=text/plain' > ${condition}_down_results.json
    else
        echo "No down regulated miRNA target genes to perform analysis on."
    fi


    """
}