#!/usr/bin/env python3

import pandas as pd
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
        "--clipExpNum",
        type=int,
        required=True,
        help="Minimum number of supporting CLIP-seq experiments required for miRNA-mRNA interaction.",
    )

    arg_parser.add_argument(
        "--degraExpNum",
        type=int,
        required=True,
        help="Minimum number of supporting degradome-seq experiments.",
    )

    arg_parser.add_argument(
        "--programNum",
        type=int,
        required=True,
        help="Minimum number of target-predicting programs. <= .7",
    )

    arg_parser.add_argument(
        "--program",
        type=str,
        required=True,
        help="Target-predicting programs (PITA, RNA22, miRmap, DIANA-microT, miRanda, PicTar, and TargetScan).",
    )

    args_parsed = arg_parser.parse_args()
    miRNA_list_file = args_parsed.miRNA_list
    deseq2_output = args_parsed.deseq2_output
    clipExpNum = args_parsed.clipExpNum
    degraExpNum = args_parsed.degraExpNum
    programNum = args_parsed.programNum
    program = args_parsed.program

    logger.info(f"miRNA_list_file: {miRNA_list_file}")
    logger.info(f"deseq2_output: {deseq2_output}")
    logger.info(f"clipExpNum: {clipExpNum}")
    logger.info(f"degraExpNum: {degraExpNum}")
    logger.info(f"programNum: {programNum}")
    logger.info(f"program: {program}")

    if deseq2_output == 1:
        # The columns are expected to follow the format:
        #   <miRNA>,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,threshold

        # Separate the miRNA_list file into 2 files: up and down regulated miRNAs
        miRNA_DE = pd.read_csv(miRNA_list_file, header=0)
        print(miRNA_DE.columns)
        print(miRNA_DE['threshold'])
        up_regulated = miRNA_DE[(miRNA_DE['log2FoldChange'] > 0) & (miRNA_DE['threshold'] == True)]
        up_regulated.columns = ['miRNA'] + up_regulated.columns[1:].tolist()
        down_regulated = miRNA_DE[(miRNA_DE['log2FoldChange'] < 0) & (miRNA_DE['threshold'] == True)]
        down_regulated.columns = ['miRNA'] + down_regulated.columns[1:].tolist()

        print(up_regulated.columns)

        # Write the up and down regulated miRNAs to a file
        up_regulated.to_csv("up_regulated_miRNAs.tsv", sep=',', index=False)
        down_regulated.to_csv("down_regulated_miRNAs.tsv", sep=',', index=False)
        
    elif deseq2_output == 0:
        print("Non DESeq2 output file being processed")
        # Header expected to follow the format: miRNA   <up or down regulated>regulated
        miRNA_list = pd.read_csv(miRNA_list_file, sep='\t', header=0)
        up_regulated = miRNA_list[miRNA_list['regulated'] == 'up']
        down_regulated = miRNA_list[miRNA_list['regulated'] == 'down']
    else:
        print("Please provide a valid value (yes or no) for the deseq2_output argument.")

    up_regulated_miRNA_targets = pd.DataFrame()
    down_regulated_miRNA_targets = pd.DataFrame()

    # Get the miRNA targets df from the ENCORI database
    up_list = []
    for index,row in up_regulated.iterrows():
        miRNA_targets = retrieve_miRNA_csvs(row['miRNA'], clipExpNum, degraExpNum, programNum, program)

        up_list.append(miRNA_targets)
    up_regulated_miRNA_targets = pd.concat(up_list)

    down_list = []
    for index,row in down_regulated.iterrows():
        miRNA_targets = retrieve_miRNA_csvs(row['miRNA'], clipExpNum, degraExpNum, programNum, program)

        down_list.append(miRNA_targets)
        
    down_regulated_miRNA_targets = pd.concat(down_list)

    # Write all the information to a csv file
    up_regulated_miRNA_targets.to_csv("all_up_regulated_ENCORI_miRNA_targets.tsv", sep='\t', index=False)
    down_regulated_miRNA_targets.to_csv("all_down_regulated_ENCORI_miRNA_targets.tsv", sep='\t', index=False)



if __name__ == "__main__":
    main()