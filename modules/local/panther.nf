process PANTHER {
    label "process_low"

    input:
    tuple file(targets_up), file(targets_down)

    output:
    tuple path("up_*_results.json"), path("down_*_results.json"), emit: panther_results

    script:    
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
    // TODO allow inputting multiple annotDataSets
    def annotDataSet = params.panther_annotDataSet ?: "GO:0008150"

    // correction options:
    // FDR - False Discovery Rate
    // BONFERRONI - Bonferroni
    // None - No correction
    def correction = params.panther_correction ?: "FDR"

    // A maximum of 100,000 genes can be submitted at a time
    """
    #!/bin/bash

    # Convert user list into array in bash
    user_input=(${annotDataSet})

    # Do input checking
    valid_annotations=("GO:0008150" "GO:0005575" "GO:0003674" \\
              "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF" "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP" \\
              "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC" "ANNOT_TYPE_ID_PANTHER_PC" \\
              "ANNOT_TYPE_ID_PANTHER_PATHWAY" "ANNOT_TYPE_ID_REACTOME_PATHWAY")
              

    valid_user_list=()
    for item in "\${user_input[@]}"; do
        if [[ " \${valid_annotations[@]} " == *" \${item} "* ]]; then
            valid_user_list+=("\${item}")
        fi
    done

    # If the valid_user_list is empty, exit the program
    if [ \${#valid_user_list[@]} -eq 0 ]; then
        echo "Invalid annotation data set. Exiting program"
        exit 1
    fi

    # Print the valid annotations that will be used
    echo "Valid annotation data sets: \${valid_user_list[@]}"

    if [[ "${correction}" == "FDR" || "${correction}" == "BONFERRONI" || "${correction}" == "None" ]]; then
        echo "Valid correction set."
    else
        echo "Invalid correction specified"
        exit 1
    fi

    # Check that the targets_up have less than 100,000 genes
    skip_up=0
    if [ \$(wc -l < "${targets_up}") -gt 100000 ]; then
        echo "Up regulated miRNA target genes file has more than 100,000 genes"
        echo "Please take the file and run directly in PANTHER website."
        skip_up=1
    fi

    # Check that the targets_down have less than 100,000 genes
    skip_down=0
    if [ \$(wc -l < ${targets_down}) -gt 100000 ]; then
        echo "Down regulated miRNA target genes file has more than 100,000 genes"
        echo "Please take the file and run directly in PANTHER website."
        skip_down=1
    fi

    echo "Input checks have been passed."

    # Iterate through the annotation datasets for up regulated genes
    for annot in "\${valid_user_list[@]}"; do
        echo "Performing analysis on up regulated genes for annotation data set: \${annot}"
        # Check if the up reglated file is empty
        if [[ -s "${targets_up}" && "\$skip_down" -eq 0 ]]; then
            curl -X POST 'https://pantherdb.org/services/oai/pantherdb/enrich/statenrich' \\
                -H 'accept: application/json' \\
                -H 'Content-Type: multipart/form-data' \\
                -F 'organism=${organism}' \\
                -F "annotDataSet=\${annot}" \\
                -F 'correction=${correction}' \\
                -F 'geneExp=@${targets_up};type=text/plain' > up_\${annot}_results.json
        else
            echo "No up regulated miRNA target genes to perform analysis on."
        fi
    done

    for annot in "\${valid_user_list[@]}"; do
        echo "Performing analysis on down regulated genes for annotation data set: \${annot}"
        # Check if the down regulated file is empty
        if [[ -s "${targets_down}" && "\$skip_down" -eq 0 ]]; then
            curl -X POST 'https://pantherdb.org/services/oai/pantherdb/enrich/statenrich' \\
                -H 'accept: application/json' \\
                -H 'Content-Type: multipart/form-data' \\
                -F 'organism=${organism}' \\
                -F "annotDataSet=\${annot}" \\
                -F 'correction=${correction}' \\
                -F 'geneExp=@${targets_down};type=text/plain' > down_\${annot}_results.json
        else
            echo "No down regulated miRNA target genes to perform analysis on."
        fi
    done

    """
}