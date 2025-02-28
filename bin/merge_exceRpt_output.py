#!/usr/bin/env python3

import logging
import pandas as pd
import argparse

logging.basicConfig(filename='merge_exceRpt_output.log', format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main():
    

    arg_parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    arg_parser.add_argument(
        "--exceRpt_results_dirs",
        type=str,
        required=True,
        help="The directory with the exceRpt results.",
    )

    arg_parsed = arg_parser.parse_args()
    exceRpt_results_dirs = arg_parsed.exceRpt_results_dirs.split(',')

    # In this script, we will combine the count results for miRNA, piRNA, tRNA, circRNA, GENCODE counts
    #       additionally, we will calculate reads per million for each sample and combine all samples

    # Create list for each type of RNA which will be converted to pandas df at the end
    miRNA_counts = []
    miRNA_rpm = []

    pirna_counts = []
    pirna_rpm = []

    trna_counts = []
    trna_rpm = []

    circrna_counts = []
    circrna_rpm = []

    gencode_counts = []
    gencode_rpm = []

    rRNA_counts = []

    # Create dataframe for a list of all of the biotypes in the samples
    biotype = pd.DataFrame(columns=['Sample_ID','miRNA','tRNA','snoRNA','circRNA','retained_intron',
                                'protein_coding','misc_RNA','Mt_tRNA','piRNA',
                                'processed_transcript','Mt_rRNA','lincRNA',
                                'antisense','snRNA','nonsense_mediated_decay',
                                'scaRNA','processed_pseudogene','sense_intronic',
                                'ribozyme','TEC','unprocessed_pseudogene',
                                'transcribed_processed_pseudogene',
                                'sense_overlapping','vaultRNA','non_coding',
                                'transcribed_unprocessed_pseudogene','macro_lncRNA',
                                'rRNA','unitary_pseudogene','non_stop_decay','sRNA',
                                '3prime_overlapping_ncrna','IG_C_gene','IG_V_gene',
                                'exogenous_genomes','TR_V_gene','pseudogene',
                                'polymorphic_pseudogene','IG_V_pseudogene',
                                'TR_J_gene','TR_C_gene','TR_V_pseudogene',
                                'IG_C_pseudogene','bidirectional_promoter_lncrna',
                                'IG_J_gene','TR_J_pseudogene','transcribed_unitary_pseudogene',
                                'translated_unprocessed_pseudogene'])

    # for each sample in the directory, read the miRNA
    for sample in exceRpt_results_dirs:
        sample_id = sample.split('/')[-1].split('_exceRpt')[0]
        logger.info(f"Processing {sample_id}")

        path_to_file = f"{sample}/{sample_id}/"

        #################
        # Mapping stats
        #################
        mapping_stats = pd.read_csv(f"{sample}/{sample_id}.stats", sep='\t', comment='-', header=None, index_col=0)
        # The ReadCount for genome is used for normalization
        genome_read_count = mapping_stats.at['genome', 1]
        logger.info(f"Genome read count: {genome_read_count}")

        #########
        # miRNA
        #########
        mirna_mature_sense = pd.read_csv(f"{path_to_file}readCounts_miRNAmature_sense.txt", sep='\t', header=0)

        # Split the ReferenceID so that the miRBase ID is the row index. Some have more than 1 miRNA in it. 
        # TODO split by | and then by :
        mirna_mature_sense['miRBase_ID'] = mirna_mature_sense['ReferenceID'].str.split(':').str[0]

        # multi-mapped adjusted counts
        temp = mirna_mature_sense[['miRBase_ID', 'multimapAdjustedReadCount']]
        temp.columns = ['miRBase_ID', sample_id]
        miRNA_counts.append(temp)

        # Get the biotype count of miRNA
        mirna_biotype = mirna_mature_sense['multimapAdjustedReadCount'].sum()

        # Calculate RPM for miRNA
        mirna_mature_sense[sample_id] = (mirna_mature_sense['multimapAdjustedReadCount'] / genome_read_count) * 1000000
        temp = mirna_mature_sense[['miRBase_ID', sample_id]]
        miRNA_rpm.append(temp)

        mirna_precursor_sense = pd.read_csv(f"{path_to_file}readCounts_miRNAprecursor_sense.txt", sep='\t', header=0)
        mirna_precursor_antisense = pd.read_csv(f"{path_to_file}readCounts_miRNAprecursor_antisense.txt", sep='\t', header=0)

        #########
        # piRNA
        #########
        pirna_sense = pd.read_csv(f"{path_to_file}readCounts_piRNA_sense.txt", sep='\t', header=0)
        pirna_antisense = pd.read_csv(f"{path_to_file}readCounts_piRNA_antisense.txt", sep='\t', header=0)

        pirna_sense['piRNA_ID'] = pirna_sense['ReferenceID'].str.split('|').str[0]

        # multi-mapped adjusted counts
        temp = pirna_sense[['piRNA_ID', 'multimapAdjustedReadCount']]
        temp.columns = ['piRNA_ID', sample_id]
        pirna_counts.append(temp)

        # Get the biotype count of piRNA
        pirna_biotype = pirna_sense["multimapAdjustedReadCount"].sum()

        # Calculate RPM for piRNA
        pirna_sense[sample_id] = (pirna_sense['multimapAdjustedReadCount'] / genome_read_count) * 1000000
        temp = pirna_sense[['piRNA_ID', sample_id]]
        pirna_rpm.append(temp)

        ########
        # tRNA
        ########
        trna_sense = pd.read_csv(f"{path_to_file}readCounts_tRNA_sense.txt", sep='\t', header=0)
        trna_antisense = pd.read_csv(f"{path_to_file}readCounts_tRNA_antisense.txt", sep='\t', header=0)

        # multi-mapped adjusted counts
        temp = trna_sense[['ReferenceID', 'multimapAdjustedReadCount']]
        temp.columns = ['tRNA_ID', sample_id]
        trna_counts.append(temp)

        # Get the biotype count of tRNA
        trna_biotype = trna_sense['multimapAdjustedReadCount'].sum()

        # Calculate RPM for tRNA
        trna_sense[sample_id] = (trna_sense['multimapAdjustedReadCount'] / genome_read_count) * 1000000
        temp = trna_sense[['ReferenceID', sample_id]]
        temp.columns = ['tRNA_ID', sample_id]
        trna_rpm.append(temp)

        ###########
        # circRNA
        ###########
        circrna_sense = pd.read_csv(f"{path_to_file}readCounts_circRNA_sense.txt", sep='\t', header=0)

        # multi-mapped adjusted counts
        temp = circrna_sense[['ReferenceID', 'multimapAdjustedReadCount']]
        temp.columns = ['circRNA_ID', sample_id]
        circrna_counts.append(temp)

        # Get the biotype count of circRNA
        circrna_biotype = circrna_sense['multimapAdjustedReadCount'].sum()

        # Calculate RPM for circRNA
        circrna_sense[sample_id] = (circrna_sense['multimapAdjustedReadCount'] / genome_read_count) * 1000000
        temp = circrna_sense[['ReferenceID', sample_id]]
        temp.columns = ['circRNA_ID', sample_id]
        circrna_rpm.append(temp)

        ##########
        # gencode
        ##########
        gencode_sense = pd.read_csv(f"{path_to_file}readCounts_gencode_sense.txt", sep='\t', header=0)
        gencode_antisense = pd.read_csv(f"{path_to_file}readCounts_gencode_antisense.txt", sep='\t', header=0)
        gencode_geneLevel_sense = pd.read_csv(f"{path_to_file}readCounts_gencode_sense_geneLevel.txt", sep='\t', header=0)
        gencode_geneLevel_antisense = pd.read_csv(f"{path_to_file}readCounts_gencode_antisense_geneLevel.txt", sep='\t', header=0)

        # multi-mapped adjusted counts
        temp = gencode_sense[['ReferenceID', 'multimapAdjustedReadCount']]
        temp.columns = ['GENCODE_ID', sample_id]
        gencode_counts.append(temp)

        # Calculate RPM for gencode
        gencode_sense[sample_id] = (gencode_sense['multimapAdjustedReadCount'] / genome_read_count) * 1000000
        temp = gencode_sense[['ReferenceID', sample_id]]
        temp.columns = ['GENCODE_ID', sample_id]
        gencode_rpm.append(temp)


        ############
        # Biotypes
        ############
        # Get the biotype counts from GENCODE
        gencode_sense['biotype'] = gencode_sense['ReferenceID'].str.split(':').str[1]
        gencode_sense_biotype = gencode_sense.groupby('biotype').sum()
        gencode_sense_biotype = gencode_sense_biotype[['multimapAdjustedReadCount']]

        # Add to miRNA biotype count
        if 'miRNA' in gencode_sense_biotype.index:
            gencode_sense_biotype.at['miRNA', 'multimapAdjustedReadCount'] += mirna_biotype
        
        gencode_sense_biotype = gencode_sense_biotype.T
        gencode_sense_biotype['Sample_ID'] = sample_id
        gencode_sense_biotype['tRNA'] = trna_biotype
        gencode_sense_biotype['circRNA'] = circrna_biotype
        gencode_sense_biotype['piRNA'] = pirna_biotype

        gencode_sense_biotype.reset_index(drop=True, inplace=True)
        biotype = pd.concat([biotype, gencode_sense_biotype]).fillna(0)

        logger.info(gencode_sense_biotype)

        #######
        # rRNA
        #######
        with open(f"{path_to_file}{sample_id}.clipped.trimmed.filtered.rRNA.readCount", 'r') as rRNA_file:
            rRNA_filtering_stats = int(rRNA_file.read().strip())
        logger.info(f"rRNA filtering stats: {rRNA_filtering_stats}")
        # Add to rRNA array
        rRNA_counts.append([sample_id, rRNA_filtering_stats])


    # Convert the lists to pandas dataframes
    merged_miRNA_counts = pd.concat(miRNA_counts, axis=0).fillna(0).groupby('miRBase_ID').sum()
    merged_miRNA_counts.to_csv('miRNA_counts.tsv', sep='\t', index=True, header=True)
    merged_miRNA_rpm = pd.concat(miRNA_rpm, axis=0).fillna(0).groupby('miRBase_ID').sum()
    merged_miRNA_rpm.to_csv('miRNA_rpm.tsv', sep='\t', index=True, header=True)

    merged_pirna_counts = pd.concat(pirna_counts, axis=0).fillna(0).groupby('piRNA_ID').sum()
    merged_pirna_counts.to_csv('piRNA_counts.tsv', sep='\t', index=True, header=True)
    merged_pirna_rpm = pd.concat(pirna_rpm, axis=0).fillna(0).groupby('piRNA_ID').sum()
    merged_pirna_rpm.to_csv('piRNA_rpm.tsv', sep='\t', index=True, header=True)

    merged_trna_counts = pd.concat(trna_counts, axis=0).fillna(0).groupby('tRNA_ID').sum()
    merged_trna_counts.to_csv('tRNA_counts.tsv', sep='\t', index=True, header=True)
    merged_trna_rpm = pd.concat(trna_rpm, axis=0).fillna(0).groupby('tRNA_ID').sum()
    merged_trna_rpm.to_csv('tRNA_rpm.tsv', sep='\t', index=True, header=True)

    merged_circrna_counts = pd.concat(circrna_counts, axis=0).fillna(0).groupby('circRNA_ID').sum()
    merged_circrna_counts.to_csv('circRNA_counts.tsv', sep='\t', index=True, header=True)
    merged_circrna_rpm = pd.concat(circrna_rpm, axis=0).fillna(0).groupby('circRNA_ID').sum()
    merged_circrna_rpm.to_csv('circRNA_rpm.tsv', sep='\t', index=True, header=True)

    merged_gencode_counts = pd.concat(gencode_counts, axis=0).fillna(0).groupby('GENCODE_ID').sum()
    merged_gencode_counts.to_csv('GENCODE_counts.tsv', sep='\t', index=True, header=True)
    merged_gencode_rpm = pd.concat(gencode_rpm, axis=0).fillna(0).groupby('GENCODE_ID').sum()
    merged_gencode_rpm.to_csv('GENCODE_rpm.tsv', sep='\t', index=True, header=True)

    biotype = biotype.T
    biotype.to_csv('biotype_counts.tsv', sep='\t', index=True, header=False)

    rRNA_df = pd.DataFrame(rRNA_counts, columns=['Sample_ID', 'rRNA_filtered_read_counts'])
    rRNA_df.to_csv('rRNA_filtered_read_counts.tsv', sep='\t', index=False, header=True)

if __name__ == "__main__":
    main()