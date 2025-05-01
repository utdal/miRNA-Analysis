#!/usr/bin/env python3

import pandas as pd
import numpy as np
import requests
import argparse
import logging

logging.basicConfig(filename='get_miRNA_targets_from_database.log', format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Function to get the data about the miRNA from ENCORI
# https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all&cellType=all'

# TargetScan has predicted miRNA targets
# ENCORI has validated miRNA targets

# Retreive from ENCORI
def retrieve_miRNA_csvs(miRNA, clipExpNum, degraExpNum, programNum, program):
    url = f'https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum={clipExpNum}&degraExpNum={degraExpNum}&pancancerNum=0&programNum={programNum}&program={program}&target=all&cellType=all'
    response = requests.get(url)
    if response.status_code == 200:
        with open(f'{miRNA}_ENCORI_targets.csv', 'wb') as f:
            f.write(response.content)
    else:
        logger.info(f'Error: {response.status_code}')
    
    # Read the csv into a dataframe and return

    miRNA_targets = pd.read_csv(f'{miRNA}_ENCORI_targets.csv', sep='\t', header=3)
    if "Or the input of" in miRNA_targets.iat[0,0]:
        logger.info(f"ENCORI does not have data on: {miRNA}")
        # Make an empty dataframe with the same columns
        miRNA_targets = pd.DataFrame(columns=miRNA_targets.columns)

    return miRNA_targets   


def main():
    arg_parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser.add_argument(
        "--miRNA_list",
        type=str,
        required=True,
        help="The file with list of differentially expressed miRNAs. Maybe the output of DESeq2 or a file with the following format: miRNA   <up or down>regulated. See documentation for more information on format.",
    )

    arg_parser.add_argument(
        "--deseq2_output",
        type=int,
        required=True,
        help="If the miRNA_list file is the output of DESeq2. yes or no.",
    )

    arg_parser.add_argument(
        "--hsa_miRTarBase_TargetScan_db",
        type=str,
        required=True,
        help="The file with the hsa_miRTarBase_TargetScan database.",
    )

    arg_parser.add_argument(
        "--experimental_evidence",
        type=str,
        required=True,
        help="How strong the experimental evidence should be. Options: 'strong', 'weak'.",
    )

    arg_parser.add_argument(
        "--min_weighted_context_percentile",
        type=float,
        required=True,
        help="The minimum weighted context percentile for predicted targets. Options: 0-100. If no predicted targets are to be used enter 101.",
    )

    arg_parser.add_argument(
        "--max_predicted_KD",
        type=float,
        required=True,
        help="The minimum predicted KD for predicted targets. Smaller KD is a stronger prediction. Options: -inf to 0. If no predicted targets are to be used enter 1.",
    )

    args_parsed = arg_parser.parse_args()
    miRNA_list_file = args_parsed.miRNA_list
    deseq2_output = args_parsed.deseq2_output
    hsa_miRTarBase_TargetScan_file = args_parsed.hsa_miRTarBase_TargetScan_db
    experimental_evidence = args_parsed.experimental_evidence
    min_weighted_context_percentile = args_parsed.min_weighted_context_percentile
    max_predicted_KD = args_parsed.max_predicted_KD

    logger.info(f"miRNA_list_file: {miRNA_list_file}")
    logger.info(f"deseq2_output: {deseq2_output}")
    logger.info(f"experimental_evidence: {experimental_evidence}")
    logger.info(f"min_weighted_context_percentile: {min_weighted_context_percentile}")
    logger.info(f"max_predicted_KD: {max_predicted_KD}")

    # Read the hsa_miRTarBase_TargetScan database
    db = pd.read_csv(hsa_miRTarBase_TargetScan_file, sep='\t', header=0)

    db['weighted context++ score percentile'].replace(".", 0, inplace=True)
    db['weighted context++ score percentile'] = db['weighted context++ score percentile'].astype(float)

    db['Predicted relative KD'].replace(".", 1, inplace=True)
    db['Predicted relative KD'] = db['Predicted relative KD'].astype(float)

    db.head(10).to_csv("db_head.tsv", sep='\t', index=False, header=True)
    logger.info(db['weighted context++ score percentile'].dtype)
    logger.info(db['Predicted relative KD'].dtype)

    if deseq2_output == 1:
        # The columns are expected to follow the format:
        #   <miRNA>,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,threshold

        # Separate the miRNA_list file into 2 files: up and down regulated miRNAs
        miRNA_DE = pd.read_csv(miRNA_list_file, header=0)
        up_regulated = miRNA_DE[(miRNA_DE['log2FoldChange'] > 0) & (miRNA_DE['threshold'] == True)]
        up_regulated.columns = ['miRNA'] + up_regulated.columns[1:].tolist()
        down_regulated = miRNA_DE[(miRNA_DE['log2FoldChange'] < 0) & (miRNA_DE['threshold'] == True)]
        down_regulated.columns = ['miRNA'] + down_regulated.columns[1:].tolist()

        # Write the up and down regulated miRNAs to a file
        up_regulated.to_csv("up_regulated_miRNAs.tsv", sep=',', index=False)
        down_regulated.to_csv("down_regulated_miRNAs.tsv", sep=',', index=False)
        
    elif deseq2_output == 0:
        logger.info("Non DESeq2 output file being processed")
        # Header expected to follow the format: miRNA   <up or down regulated>regulated
        miRNA_list = pd.read_csv(miRNA_list_file, sep='\t', header=0)
        if len(miRNA_list.columns) == 1:
            up_regulated = miRNA_list
            down_regulated = pd.DataFrame(columns=miRNA_list.columns)
        else:
            up_regulated = miRNA_list[miRNA_list['regulated'] == 'up']
            down_regulated = miRNA_list[miRNA_list['regulated'] == 'down']
    else:
        logger.info("Please provide a valid value (yes or no) for the deseq2_output argument.")

    up_regulated_miRNA_targets = pd.DataFrame()
    down_regulated_miRNA_targets = pd.DataFrame()

    # Get the miRNA targets df from the database
    up_mirna = up_regulated['miRNA'].unique()
    down_mirna = down_regulated['miRNA'].unique()

    up_regulated_miRNA_targets = db[db['miRNA'].isin(up_mirna)]
    down_regulated_miRNA_targets = db[db['miRNA'].isin(down_mirna)]

    # Set the Kd column to a numeric value
    up_regulated_miRNA_targets['Predicted relative KD'] = pd.to_numeric(up_regulated_miRNA_targets['Predicted relative KD'], errors='coerce')
    up_regulated_miRNA_targets['Predicted relative KD'] = up_regulated_miRNA_targets['Predicted relative KD'].fillna(float('inf'))
    down_regulated_miRNA_targets['Predicted relative KD'] = pd.to_numeric(down_regulated_miRNA_targets['Predicted relative KD'], errors='coerce')
    down_regulated_miRNA_targets['Predicted relative KD'] = down_regulated_miRNA_targets['Predicted relative KD'].fillna(float('inf'))

    if experimental_evidence == 'strong':
        up_regulated_miRNA_targets = up_regulated_miRNA_targets[(up_regulated_miRNA_targets['Experimental Evidence'] == experimental_evidence) | 
                                                                (up_regulated_miRNA_targets['weighted context++ score percentile'] >= min_weighted_context_percentile) & 
                                                                (up_regulated_miRNA_targets['Predicted relative KD'] <= max_predicted_KD)]
        down_regulated_miRNA_targets = down_regulated_miRNA_targets[(down_regulated_miRNA_targets['Experimental Evidence'] == experimental_evidence) | 
                                                                    (down_regulated_miRNA_targets['weighted context++ score percentile'] >= min_weighted_context_percentile) & 
                                                                    (down_regulated_miRNA_targets['Predicted relative KD'] <= max_predicted_KD)]
    else:
        up_regulated_miRNA_targets = up_regulated_miRNA_targets[(up_regulated_miRNA_targets['Experimental Evidence'] == 'weak') | 
                                                                up_regulated_miRNA_targets['Experimental Evidence'] == 'strong' | 
                                                                (up_regulated_miRNA_targets['weighted context++ score percentile'] >= min_weighted_context_percentile) & 
                                                                (up_regulated_miRNA_targets['Predicted relative KD'] <= max_predicted_KD)]
        down_regulated_miRNA_targets = down_regulated_miRNA_targets[(down_regulated_miRNA_targets['Experimental Evidence'] == 'weak') | 
                                                                    down_regulated_miRNA_targets['Experimental Evidence'] == 'strong' | 
                                                                    (down_regulated_miRNA_targets['weighted context++ score percentile'] >= min_weighted_context_percentile) & 
                                                                    (down_regulated_miRNA_targets['Predicted relative KD'] <= max_predicted_KD)]

    # Convert the inf back to .
    up_regulated_miRNA_targets['Predicted relative KD'] = up_regulated_miRNA_targets['Predicted relative KD'].replace(float('inf'), '.')
    down_regulated_miRNA_targets['Predicted relative KD'] = down_regulated_miRNA_targets['Predicted relative KD'].replace(float('inf'), '.')

    # Write all the information to a csv file
    up_regulated_miRNA_targets.to_csv("all_up_reg_miRNA_targets_miRTarBase_TargetScan.tsv", sep='\t', index=False)
    down_regulated_miRNA_targets.to_csv("all_down_reg_miRNA_targets_miRTarBase_TargetScan.tsv", sep='\t', index=False)

    logger.info("Done")

if __name__ == "__main__":
    main()