import pandas as pd
import requests
import argparse

# Function to get the data about the miRNA from ENCORI
# https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all&cellType=all'

# TargetScan has predicted miRNA targets
# ENCORI has validated miRNA targets

def retrieve_miRNA_csvs(miRNA, regulated, clipExpNum, degraExpNum, programNum, program):
    url = f'https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum={clipExpNum}&degraExpNum={degraExpNum}&pancancerNum=0&programNum={programNum}&program={program}&target=all&cellType=all'
    response = requests.get(url)
    if response.status_code == 200:
        with open(f'{miRNA}_ENCORI_targets.csv', 'wb') as f:
            f.write(response.content)
    else:
        print(f'Error: {response.status_code}')
    
    # Read the csv into a dataframe and return
    miRNA_targets = pd.read_csv(f'{miRNA}_ENCORI_targets.csv', sep='\t', header=3)
    miRNA_targets['regulation'] = regulated

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
        help="The file with list of differentially expressed miRNAs.",
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
    miRNA_list = args_parsed.miRNA_list
    clipExpNum = args_parsed.clipExpNum
    degraExpNum = args_parsed.degraExpNum
    programNum = args_parsed.programNum
    program = args_parsed.program

    # Make a dataframe for the all targets of the differentially expressed miRNAs
    all_validated_miRNA_targets = pd.DataFrame()

    # Header expected to follow the format: miRNA   <up or down regulated>regulated
    miRNA_list = pd.read_csv(miRNA_list, sep='\t', header=0)

    # Get the miRNA targets df from the ENCORI database
    for index,row in miRNA_list.iterrows():
        miRNA_targets = retrieve_miRNA_csvs(row['miRNA'], row['regulated'], clipExpNum, degraExpNum, programNum, program)

        # TODO Try to reduce complexity
        temp = all_validated_miRNA_targets
        all_validated_miRNA_targets = pd.concat([temp,miRNA_targets])

    # Write all the information to a csv file
    all_validated_miRNA_targets.to_csv("all_ENCORI_miRNA_DE_targets.csv", sep='\t')



if __name__ == "__main__":
    main()