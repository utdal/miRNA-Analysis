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
        required=True,
        help="The file with the bulk RNA-seq counts. Format should be gene_name,sample1,sample2,...,sampleN in header.",
    )

    arg_parser.add_argument(
        "--min_expression",
        type=int,
        required=False,
        default=50,
        help="The minimum expression of the gene to be considered.",
    )

    arg_parser.add_argument(
        "--condition",
        type=str,
        required=True,
        help="The condition differential expressed miRNA.",
    )

    # Check paths are correct
    args_parsed = arg_parser.parse_args()
    up_miRNA_targets_file = args_parsed.up_miRNA_targets
    down_miRNA_targets_file = args_parsed.down_miRNA_targets
    bulkRNAseq_file = args_parsed.bulkRNAseq
    min_expression = args_parsed.min_expression
    condition = args_parsed.condition

    up_miRNA_targets = pd.read_csv(up_miRNA_targets_file, sep='\t', header=0)
    down_miRNA_targets = pd.read_csv(down_miRNA_targets_file, sep='\t', header=0)
    bulkRNAseq = pd.read_csv(bulkRNAseq_file, sep=',', header=0)

    # Header of the miRNA_targets file
    # miRNAid	miRNAname	geneID	geneName	geneType	chromosome	narrowStart	narrowEnd	broadStart	broadEnd	strand	clipExpNum	degraExpNum	RBP	PITA	RNA22	miRmap	microT	miRanda	PicTar	TargetScan	TDMDScore	phyloP	pancancerNum	cellline/tissue	regulation

    # Get the average gene counts across all samples
    bulkRNAseq['mean'] = bulkRNAseq.iloc[:,1:].mean(axis=1)

    # Filter the genes with a mean expression higher than the minimum expression
    bulkRNAseq_min = bulkRNAseq[bulkRNAseq['mean'] > min_expression]

    logger.info("bulkRNAseq file head with mean:")
    logger.info(len(bulkRNAseq_min))

    # Filter the up regulated miRNA targets
    filtered_up_miRNA_targets = up_miRNA_targets[up_miRNA_targets['geneID'].isin(bulkRNAseq_min['gene_name']) | up_miRNA_targets["geneName"].isin(bulkRNAseq_min['gene_name'])]
    filtered_up_miRNA_targets.to_csv(f"{condition}_filtered_ENCORI_up_miRNA_DE_targets.tsv", sep='\t', index=False)

    # Filter the down regulated miRNA targets
    filtered_down_miRNA_targets = down_miRNA_targets[down_miRNA_targets['geneID'].isin(bulkRNAseq_min['gene_name']) | down_miRNA_targets["geneName"].isin(bulkRNAseq_min['gene_name'])]
    filtered_down_miRNA_targets.to_csv(f"{condition}_filtered_ENCORI_down_miRNA_DE_targets.tsv", sep='\t', index=False)

    # Get the number of miRNAs targeting each gene within the filtered tissue gene expression
    up_regulated_final_target_gene_targeted_count = (
        filtered_up_miRNA_targets.groupby('geneName')['miRNAname']
        .nunique()
        .reset_index()
        .rename(columns={'miRNAname': 'miRNA_DE_targets_count'})
    )

    down_regulated_final_target_gene_targeted_count = (
        filtered_down_miRNA_targets.groupby('geneName')['miRNAname']
        .nunique()
        .reset_index()
        .rename(columns={'miRNAname': 'miRNA_DE_targets_count'})
    )

    up_regulated_final_target_gene_targeted_count.to_csv(f"{condition}_up_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)
    down_regulated_final_target_gene_targeted_count.to_csv(f"{condition}_down_regulated_final_target_gene_targeted_count.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()