#!/usr/bin/env python3

import json
import pandas as pd
import logging
import argparse

# For graphing
import os
os.environ['MPLCONFIGDIR'] = '/tmp'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap


def main():
    print("Starting main function")
    logging.basicConfig(filename='interpreting_panther.log', format='%(asctime)s %(message)s', filemode='w')
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    arg_parser = argparse.ArgumentParser(
        description="Script to convert PANTHER output into a readable format.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser.add_argument(
        "--up_panther_output",
        type=str,
        required=True,
        help="The file with the PANTHER output of up regulated miRNAs' targets.",
    )

    arg_parser.add_argument(
        "--down_panther_output",
        type=str,
        required=True,
        help="The file with the PANTHER output of down regulated miRNAs' targets.",
    )

    args_parsed = arg_parser.parse_args()
    up_panther_output_files = args_parsed.up_panther_output.split(" ")
    down_panther_output_files = args_parsed.down_panther_output.split(" ")

    for up_panther_output_file, down_panther_output_file in zip(up_panther_output_files, down_panther_output_files):
        logger.info(f"Processing up panther output file: {up_panther_output_file} and down panther output file: {down_panther_output_file}")
        # Annotation Data
        annot = up_panther_output_file.split("up_")[1].split("_results")[0]

        os.mkdir(annot)

        # Get up panther file
        up_panther = pd.read_json(up_panther_output_file)
        up_panther = up_panther.T

        # Get results row into readable table format
        up_panther_results = pd.json_normalize(up_panther['result'].explode())
        up_panther_results.columns = ['Number of Genes', 'FDR', 'p-value', 'plus_minus', 'term.id', 'term.label']
        up_panther_results['GO Term'] = up_panther_results['term.label'] + " (" + up_panther_results['term.id'] + ")"
        up_panther_results = up_panther_results.sort_values(by=['FDR', 'p-value','Number of Genes'], ascending=[True, True, False])
        up_panther_results['p-value'] = up_panther_results['p-value'].astype(float)
        up_panther_results['p-value'] = up_panther_results['p-value'].apply(lambda x: np.format_float_scientific(x, precision=4))
        up_panther_results['FDR'] = up_panther_results['FDR'].astype(float)
        up_panther_results['FDR'] = up_panther_results['FDR'].apply(lambda x: np.format_float_scientific(x, precision=4))
        up_panther_results.to_csv(f"{annot}/up_enrichment_results.csv", sep="\t", header=True, index=False)

        # Get down panther file
        down_panther = pd.read_json(down_panther_output_file)
        down_panther = down_panther.T

        # Get results row into readable table format
        down_panther_results = pd.json_normalize(down_panther['result'].explode())
        down_panther_results.columns = ['Number of Genes', 'FDR', 'p-value', 'plus_minus', 'term.id', 'term.label']
        down_panther_results['GO Term'] = down_panther_results['term.label'] + " (" + down_panther_results['term.id'] + ")"
        down_panther_results = down_panther_results.sort_values(by=['FDR', 'p-value', 'Number of Genes'], ascending=[True, True, False])
        down_panther_results['p-value'] = down_panther_results['p-value'].astype(float)
        down_panther_results['p-value'] = down_panther_results['p-value'].apply(lambda x: np.format_float_scientific(x, precision=4))
        down_panther_results['FDR'] = down_panther_results['FDR'].astype(float)
        down_panther_results['FDR'] = down_panther_results['FDR'].apply(lambda x: np.format_float_scientific(x, precision=4))
        down_panther_results.to_csv(f"{annot}/down_enrichment_results.csv", sep="\t", header=True, index=False)

        # TODO split the terms in the results into clusters with similar functions
        # Sort by smalled fdr and p-value
        # Plot the ones that have an fdr < 0.1 and p-value < 0.1
        # Only perform the clustering if there are more than 30 terms

        up_clusters = split_into_clusters(up_panther_results[(up_panther_results['FDR'].astype(float) < 0.05) & (up_panther_results['p-value'].astype(float) < 0.05)])
        down_clusters = split_into_clusters(down_panther_results[(down_panther_results['FDR'].astype(float) < 0.05) & (down_panther_results['p-value'].astype(float) < 0.05)])

        # Set basic plot settings
        label_size = 12
        mpl.rcParams['xtick.labelsize'] = label_size
        mpl.rcParams['ytick.labelsize'] = label_size
        mpl.rcParams["font.sans-serif"] = ["Arial"]
        mpl.rcParams["font.family"] = "Arial"
        mpl.rcParams['pdf.fonttype'] = 42
        sns.set_theme(style="white")
        sns.set_theme(style="ticks")
        sns.set_style("whitegrid")

        # Plot each cluster for up
        cluster_count = 1
        for cluster in up_clusters:
            logger.info(cluster.describe())
            up_graph = sns.relplot(x="Number of Genes", y="GO Term", hue="FDR", 
                                size='p-value',sizes=(500, 100), legend='full', 
                                palette="Wistia", height=10, aspect=1.25, data=cluster)
            
            # Wrap the y-axis labels
            wrapped_labels = [textwrap.fill(str(label), width=30) for label in cluster['GO Term']]
            up_graph.ax.set_yticklabels(wrapped_labels)

            up_graph.savefig(f"{annot}/up_cluster{cluster_count}.png")
            logger.info(f'Finished processing up cluster {cluster_count}')
            cluster_count += 1

        # Plot each cluster for down
        cluster_count = 1
        for cluster in down_clusters:
            logger.info(cluster.describe())
            down_graph = sns.relplot(x="Number of Genes", y="GO Term", hue="FDR", 
                                size='p-value',sizes=(500, 100), legend='full', 
                                palette="Wistia", height=10, aspect=1.25, data=cluster)
            
            # Wrap the y-axis labels
            wrapped_labels = [textwrap.fill(str(label), width=30) for label in cluster['GO Term']]
            up_graph.ax.set_yticklabels(wrapped_labels)

            down_graph.savefig(f"{annot}/down_cluster{cluster_count}.png")
            logger.info(f'Finished processing down cluster {cluster_count}')
            cluster_count += 1

        logger.info("Finished processing files")
        # TODO Take rest of up info in panther into a readable format to output
        # TODO Take rest of down info in panther into a readable format to output

def split_into_clusters(df):
    # Split the terms into clusters
    return [df.iloc[i:i + 15] for i in range(0, len(df), 15)]

if __name__ == "__main__":
    main()