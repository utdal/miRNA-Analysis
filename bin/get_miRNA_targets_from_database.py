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
        type=str,
        required=True,
        help="If the miRNA_list file is the output of DESeq2.",
    )

    arg_parser.add_argument(
        "--condition",
        type=str,
        required=True,
        help="The condition from the meta2.condition.",
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
    condition = args_parsed.condition
    clipExpNum = args_parsed.clipExpNum
    degraExpNum = args_parsed.degraExpNum
    programNum = args_parsed.programNum
    program = args_parsed.program

    logger.info(f"miRNA_list: {miRNA_list}")

    # Make a dataframe for the all targets of the differentially expressed miRNAs
    all_validated_miRNA_targets = pd.DataFrame()

    if deseq2_output == 'True':
        # The miRNA_list file should be the output file from DESeq2
        #   Therefore, the columns be 
        #   <miRNA>,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj

        # Separate the miRNA_list file into 2 files: up and down regulated miRNAs
        miRNA_DE = pd.read_csv(miRNA_list_file, sep='\t', header=0)
        up_regulated = miRNA_DE[miRNA_DE['log2FoldChange'] > 0]
        up_regulated.columns = ['miRNA'] + up_regulated.columns[1:].tolist()
        down_regulated = miRNA_DE[miRNA_DE['log2FoldChange'] < 0]
        down_regulated.columns = ['miRNA'] + down_regulated.columns[1:].tolist()

        # Write the up and down regulated miRNAs to a file
        up_regulated.to_csv(f"{condition}_up_regulated_miRNAs.tsv", sep=',', index=False)
        down_regulated.to_csv(f"{condition}_down_regulated_miRNAs.tsv", sep=',', index=False)
        
    elif deseq2_output == 'False':
        # Header expected to follow the format: miRNA   <up or down regulated>regulated
        miRNA_list = pd.read_csv(miRNA_list, sep='\t', header=0)
        up_regulated = miRNA_list[miRNA_list['regulated'] == 'up']
        down_regulated = miRNA_list[miRNA_list['regulated'] == 'down']
    else:
        logger.info("Please provide a valid value (True or False) for the deseq2_output argument.")

    # Get the miRNA targets df from the ENCORI database
    for index,row in up_regulated.iterrows():
        miRNA_targets = retrieve_miRNA_csvs(row['miRNA'], clipExpNum, degraExpNum, programNum, program)

        # TODO Try to reduce complexity
        temp = up_regulated_miRNA_targets
        up_regulated_miRNA_targets = pd.concat([temp,miRNA_targets])

    for index,row in down_regulated.iterrows():
        miRNA_targets = retrieve_miRNA_csvs(row['miRNA'], clipExpNum, degraExpNum, programNum, program)

        # TODO Try to reduce complexity
        temp = down_regulated_miRNA_targets
        down_regulated_miRNA_targets = pd.concat([temp,miRNA_targets])

    # Write all the information to a csv file
    up_regulated_miRNA_targets.to_csv(f"{condition}_all_up_regulated_ENCORI_miRNA_targets.tsv", sep='\t', index=False)
    down_regulated_miRNA_targets.to_csv(f"{condition}_all_down_regulated_ENCORI_miRNA_targets.tsv", sep='\t', index=False)



if __name__ == "__main__":
    main()