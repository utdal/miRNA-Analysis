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
    all_samples.to_csv('all_samples_novel_miRNAs.csv', sep='\t', index=False)

    # Set to keep track of the visited indices
    visited_indices = set()

    # Temp arrays for memory efficiency
    similar_hairpins_array = []
    only_one_detection_array = []

    ################################################
    # SIMILAR HAIRPINS DETECTION with 90% overlap
    ################################################
    # Iterate over the rows of the dataframe to find similar hairpins
    for test_index,test_row in all_samples.iterrows():
        if test_index in visited_indices:
            continue

        has_sim = False
        for index,row in all_samples.iterrows():
            if index != test_index and index not in visited_indices and sim_hairpin_overlap(test_row['consensus precursor sequence'], row['consensus precursor sequence']):
                # copy row into similar_hairpin + delete row
                has_sim = True
                similar_hairpins_array.append(row)
                visited_indices.add(index)

        if has_sim is True:
            # write test_row to the similarl_hairpin
            similar_hairpins_array.append(test_row)
            has_sim = False
        else:
            # write test_row to only_one
            only_one_detection_array.append(test_row)

    similar_hairpins = pd.DataFrame(similar_hairpins_array, columns=all_samples.columns)
    only_one_detection = pd.DataFrame(only_one_detection_array, columns=all_samples.columns)

    similar_hairpins.to_csv('similarSort_hairpin_novel_miRNA.csv', sep='\t', index=False)
    only_one_detection.to_csv('similarSort_only_one_novel_miRNA_detected.csv', sep='\t', index=False)

    #########################################
    # EXACT HAIRPIN AND MATURE SEQ MATCHES
    #########################################

    # Reset variables
    similar_hairpins_array = []
    only_one_detection_array = []
    visited_indices = set()
    # Iterate over the similar hairpins to find exact matches
    for test_index,test_row in all_samples.iterrows():
        if test_index in visited_indices:
            continue

        has_sim = False
        for index,row in all_samples.iterrows():
            if index != test_index and index not in visited_indices and same_hairpin_mature(test_row['consensus precursor sequence'], row['consensus precursor sequence'], test_row['consensus mature sequence'], row['consensus mature sequence']):
                # copy row into similar_hairpin + delete row
                has_sim = True
                similar_hairpins_array.append(row)
                visited_indices.add(index)

        if has_sim is True:
            # write test_row to the similarl_hairpin
            similar_hairpins_array.append(test_row)
            has_sim = False
        else:
            # write test_row to only_one
            only_one_detection_array.append(test_row)

    similar_hairpins = pd.DataFrame(similar_hairpins_array, columns=all_samples.columns)
    only_one_detection = pd.DataFrame(only_one_detection_array, columns=all_samples.columns)
    similar_hairpins.to_csv('exactSort_hairpin_mature_novel_miRNA.csv', sep='\t', index=False)
    only_one_detection.to_csv('exactSort_only_one_novel_miRNA_detected.csv', sep='\t', index=False)



# Function of exact hairpin and mature seq matches
def same_hairpin_mature(hairpin1, hairpin2, mature1, mature2):
    return hairpin1 == hairpin2 and mature1 == mature2

# Function for calculated ration of overlap between 2 seqs
def overlap_ratio(seq1, seq2):
    # Identify the longest common substring between seq1 and seq2.
    # There are efficient algorithms for this (e.g., using dynamic programming),
    # but here's a basic conceptual example.

    max_overlap = 0
    len1, len2 = len(seq1), len(seq2)
    
    # Check overlap with seq1 suffix and seq2 prefix
    for i in range(1, min(len1, len2) + 1):
        if seq1[-i:] == seq2[:i]:
            max_overlap = max(max_overlap, i)
    
    # Check overlap with seq2 suffix and seq1 prefix
    for i in range(1, min(len1, len2) + 1):
        if seq2[-i:] == seq1[:i]:
            max_overlap = max(max_overlap, i)
    
    # Return ratio based on the shorter sequence length
    return max_overlap / min(len1, len2)

# Function for returning similarity based on threshold overlap ratio
def sim_hairpin_overlap(seq1, seq2, threshold=0.9):
    return overlap_ratio(seq1, seq2) >= threshold
    
if __name__ == "__main__":
    main()