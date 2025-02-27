#!/usr/bin/env python3

import pandas as pd
import io
import argparse
import logging

logging.basicConfig(filename='miRDeep2_filtering.log', format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
    arg_parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser.add_argument(
        "--miRDeep2_results_list",
        type=str,
        required=True,
        help="The directory with the miRDeep2 results.",
    )

    arg_parser.add_argument(
        "--min_miRDeep2_score",
        type=float,
        required=False,
        default=4.0,
        help="The minimum miRDeep2 score to be considered.",
    )
    args_parsed = arg_parser.parse_args()

    mirdeep_results_list = args_parsed.miRDeep2_results_list.split(',')
    min_miRDeep2_score = args_parsed.min_miRDeep2_score

    logger.info(f"All mirdeep2 results file list: {mirdeep_results_list}")

    # List to hold the processed dataframe with only the novel miRNA detections for each sample
    novel_results_samples = []
    # Iterate over the list of samples and process the miRDeep2 results in to usable format
    for file in mirdeep_results_list:
        # Filter the list for only those that have the csv file extension
        if not file.endswith(".csv"):
            continue
        
        sample_id = file[7:-4]
        logger.info(f"Processing {sample_id}: {file}")
        
        # Get only the second table from the csv file
        with open(file, 'r') as file:
            sections = file.read().strip().split('\n\n\n')
        novel_results = pd.read_csv(io.StringIO(sections[1]), sep='\t', header=1)

        # Process the dataframe
        novel_results.insert(0, "Sample Name", sample_id)

        # Only further use data with a miRDeep2 score of at least the minimum mirdeep2 score
        novel_results = novel_results.drop(novel_results[novel_results['miRDeep2 score'] <= min_miRDeep2_score].index)

        # Drop columns that are not needed
        novel_results.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)

        if novel_results.empty:
            logger.info(f"No novel miRNA detected in {sample_id}")

        # Append the processed dataframe to the list
        novel_results_samples.append(novel_results)

    # Concatenate all the dataframes in the list
    all_samples = pd.concat(novel_results_samples)
    all_samples.to_csv('all_samples.csv', sep='\t')

    # Create dataframe for the detections with similar hairpins and the ones with only one detection
    similar_hairpins = pd.DataFrame(columns=all_samples.columns)
    only_one_detection = pd.DataFrame(columns=all_samples.columns)

    # Set to keep track of the visited indices
    visited_indices = set()

    # Iterate over the rows of the all_samples dataframe
    for test_index,test_row in all_samples.iterrows():
        if test_index in visited_indices:
            continue

        has_sim = False
        for index,row in all_samples.iterrows():
            if index != test_index and index not in visited_indices and sim_hairpin_test(test_row['consensus precursor sequence'], row['consensus precursor sequence']):
                # copy row into similar_hairpin + delete row
                has_sim = True
                similar_hairpins.loc[len(similar_hairpins)] = row
                visited_indices.add(index)

        if has_sim is True:
            # write test_row to the similarl_hairpin
            similar_hairpins.loc[len(similar_hairpins)] = test_row
            has_sim = False
        else:
            # write test_row to only_one
            only_one_detection.loc[len(only_one_detection)] = test_row

    similar_hairpins.to_csv('similar_hairpin_all_samples.csv', sep='\t')

    only_one_detection.to_csv('only_one_novel_miRNA_detected_across_all_samples.csv', sep='\t')


def sim_hairpin_test(hairpin1, hairpin2):
    if hairpin1 in hairpin2 or hairpin2 in hairpin1:
            return True
    else:
        return False
    
if __name__ == "__main__":
    main()