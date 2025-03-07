#!/usr/bin/env python3

import pandas as pd
import argparse
import logging

def main():

    logging.basicConfig(filename='intersect_targets.log', format='%(asctime)s %(message)s', filemode='w')
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    arg_parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    
    arg_parser.add_argument(
        "--up_miRNA_targets",
        type=str,
        required=True,
        help="The file with the targets of up regulated miRNAs.",
    )
    arg_parser.add_argument(
        "--down_miRNA_targets",
        type=str,
        required=True,
        help="The file with the targets of down regulated miRNAs.",
    )

    arg_parser.add_argument(
        "--bulkRNAseq",
        type=str,
        required=False,
        help="The file with the bulk RNA-seq counts. Format should be gene_name,sample1,sample2,...,sampleN in header.",
    )

    arg_parser.add_argument(
        "--min_expression",
        type=int,
        required=False,
        default=50,
        help="The minimum expression of the gene to be considered.",
    )

    # Check paths are correct
    args_parsed = arg_parser.parse_args()
    up_miRNA_targets_file = args_parsed.up_miRNA_targets
    down_miRNA_targets_file = args_parsed.down_miRNA_targets
    bulkRNAseq_file = args_parsed.bulkRNAseq
    min_expression = args_parsed.min_expression

    # Read miRNA target files
    up_miRNA_targets = pd.read_csv(up_miRNA_targets_file, sep='\t', header=0)
    up_miRNA_targets['Target Gene (Entrez ID)'] = up_miRNA_targets['Target Gene (Entrez ID)'].astype(int)
    down_miRNA_targets = pd.read_csv(down_miRNA_targets_file, sep='\t', header=0)
    down_miRNA_targets['Target Gene (Entrez ID)'] = down_miRNA_targets['Target Gene (Entrez ID)'].astype(int)

    print(bulkRNAseq_file)

    # Header of the miRNA_targets file
    # miRNAid	miRNAname	geneID	geneName	geneType	chromosome	narrowStart	narrowEnd	broadStart	broadEnd	strand	clipExpNum	degraExpNum	RBP	PITA	RNA22	miRmap	microT	miRanda	PicTar	TargetScan	TDMDScore	phyloP	pancancerNum	cellline/tissue	regulation

    ################################################################################################
    # Check is if bulk RNA-seq file is provided, filter the miRNA targets by the gene expression
    ################################################################################################
    if bulkRNAseq_file == "NO_FILE_PROVIDED":
        logger.info("No bulk RNA-seq file provided, only miRNA targets will be filtered.")
    else:
        bulkRNAseq = pd.read_csv(bulkRNAseq_file, sep=',', header=0)
        # Get the average gene counts across all samples
        bulkRNAseq['mean'] = bulkRNAseq.iloc[:,1:].mean(axis=1)
        # Filter the genes with a mean expression higher than the minimum expression
        bulkRNAseq_min = bulkRNAseq[bulkRNAseq['mean'] > min_expression]

        logger.info("bulkRNAseq file head with mean:")
        logger.info(len(bulkRNAseq_min))

        # Filter the up regulated miRNA targets
        up_miRNA_targets = up_miRNA_targets[up_miRNA_targets['Target Gene (Entrez ID)'].isin(bulkRNAseq_min['gene_name']) | up_miRNA_targets["Target Gene Symbol"].isin(bulkRNAseq_min['gene_name'])]
        up_miRNA_targets.to_csv("filtered_up_miRNA_DE_targets.tsv", sep='\t', index=False)

        # Filter the down regulated miRNA targets
        down_miRNA_targets = down_miRNA_targets[down_miRNA_targets['Target Gene (Entrez ID)'].isin(bulkRNAseq_min['gene_name']) | down_miRNA_targets["Target Gene Symbol"].isin(bulkRNAseq_min['gene_name'])]
        down_miRNA_targets.to_csv("filtered_down_miRNA_DE_targets.tsv", sep='\t', index=False)

    #############################################################################################
    # Get the number of miRNAs targeting each gene within the filtered tissue gene expression
    #############################################################################################
    up_regulated_final_target_gene_targeted_count = (
        up_miRNA_targets.groupby('Target Gene Symbol')['miRNA']
        .nunique()
        .reset_index()
        .rename(columns={'miRNA': 'miRNA_DE_targets_count'})
    )

    logger.info(f"Up Regulated: Total number of interactions = {up_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'].sum()}")

    down_regulated_final_target_gene_targeted_count = (
        down_miRNA_targets.groupby('Target Gene Symbol')['miRNA']
        .nunique()
        .reset_index()
        .rename(columns={'miRNA': 'miRNA_DE_targets_count'})
    )

    logger.info(f"Down Regulated: Total number of interactions = {down_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'].sum()}")

    up_regulated_final_target_gene_targeted_count.to_csv("up_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)
    down_regulated_final_target_gene_targeted_count.to_csv("down_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()