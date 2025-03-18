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
        help="The file with the bulk RNA-seq counts. Header format must include Gene_Name anywhere.",
    )

    arg_parser.add_argument(
        "--min_expression",
        type=int,
        required=False,
        default=0,
        help="The minimum expression of the gene to be considered.",
    )

    arg_parser.add_argument(
        "--min_miRNA_target_mRNA",
        type=int,
        required=False,
        default=5,
        help="The minimum number of miRNA targeting a mRNA to be considered.",
    )

    # Check paths are correct
    args_parsed = arg_parser.parse_args()
    up_miRNA_targets_file = args_parsed.up_miRNA_targets
    down_miRNA_targets_file = args_parsed.down_miRNA_targets
    bulkRNAseq_file = args_parsed.bulkRNAseq
    min_expression = args_parsed.min_expression
    min_miRNA_target_mRNA = args_parsed.min_miRNA_target_mRNA

    # Read miRNA target files
    up_miRNA_targets = pd.read_csv(up_miRNA_targets_file, sep='\t', header=0, dtype=str)
    down_miRNA_targets = pd.read_csv(down_miRNA_targets_file, sep='\t', header=0, dtype=str)

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

        # TODO - Check if the header of the bulk RNA-seq file has 'Gene_Name' in it and counts

        # Get the average gene counts across all samples
        # Select only the numerical columns
        counts = bulkRNAseq.select_dtypes(include=['int64', 'float64'])
        bulkRNAseq['mean'] = counts.mean(axis=1)
        # Filter the genes with a mean expression higher than the minimum expression
        bulkRNAseq_min = bulkRNAseq[bulkRNAseq['mean'] > min_expression]

        logger.info("bulkRNAseq file head with mean:")
        logger.info(len(bulkRNAseq_min))

        # Filter the up regulated miRNA targets
        up_miRNA_targets = up_miRNA_targets[up_miRNA_targets['Target Gene (Entrez ID)'].isin(bulkRNAseq_min['Gene_Name']) | up_miRNA_targets["Target Gene Symbol"].isin(bulkRNAseq_min['Gene_Name'])]
        up_miRNA_targets.to_csv("RNAseq_filtered_up_miRNA_DE_targets.tsv", sep='\t', index=False)

        # Filter the down regulated miRNA targets
        down_miRNA_targets = down_miRNA_targets[down_miRNA_targets['Target Gene (Entrez ID)'].isin(bulkRNAseq_min['Gene_Name']) | down_miRNA_targets["Target Gene Symbol"].isin(bulkRNAseq_min['Gene_Name'])]
        down_miRNA_targets.to_csv("RNAseq_filtered_down_miRNA_DE_targets.tsv", sep='\t', index=False)

    #############################################################################################
    # Get the number of miRNAs targeting each gene within the filtered tissue gene expression
    #############################################################################################
    up_regulated_final_target_gene_targeted_count = (
        up_miRNA_targets.groupby('Target Gene Symbol')['miRNA']
        .agg(miRNA_DE_targets_count='nunique', miRNA_list=lambda x: '/'.join(x))
        .reset_index()
    )

    logger.info(f"Up Regulated: Total number of interactions = {up_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'].sum()}")

    down_regulated_final_target_gene_targeted_count = (
        down_miRNA_targets.groupby('Target Gene Symbol')['miRNA']
        .agg(miRNA_DE_targets_count='nunique', miRNA_list=lambda x: '/'.join(x))
        .reset_index()
    )

    logger.info(f"Down Regulated: Total number of interactions = {down_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'].sum()}")

    # Sort from highest to lowest
    up_regulated_final_target_gene_targeted_count = up_regulated_final_target_gene_targeted_count.sort_values(by='miRNA_DE_targets_count', ascending=False)
    down_regulated_final_target_gene_targeted_count = down_regulated_final_target_gene_targeted_count.sort_values(by='miRNA_DE_targets_count', ascending=False)

    up_regulated_final_target_gene_targeted_count.to_csv("up_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)
    down_regulated_final_target_gene_targeted_count.to_csv("down_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)

    if min_miRNA_target_mRNA > 0:
        up_regulated_final_target_gene_targeted_count = up_regulated_final_target_gene_targeted_count[up_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'] >= min_miRNA_target_mRNA]
        down_regulated_final_target_gene_targeted_count = down_regulated_final_target_gene_targeted_count[down_regulated_final_target_gene_targeted_count['miRNA_DE_targets_count'] >= min_miRNA_target_mRNA]

        up_regulated_final_target_gene_targeted_count.to_csv("filtered_up_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)
        down_regulated_final_target_gene_targeted_count.to_csv("filtered_down_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()