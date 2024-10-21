#!/usr/bin/env python3

print("We are here.")

from Bio import SeqIO
import re
import requests

# miRDeep2 requires that the miRNA ref files do not have letters other than AGTUCN in the sequence.
# for best human specific results, filter the for files:
#   mature homo sapiens only
#   mature not homo sapiens only
#   hairpin homo spaiens only

# Link as of 10/10/2024
# hairpin: https://www.mirbase.org/download/hairpin.fa
# mature: https://www.mirbase.org/download/mature.fa

# Link as 10/21/2024
# miRNA locations: https://www.mirbase.org/download/hsa.gff3


def download_miRBase():
    hairpin_url = 'https://www.mirbase.org/download/hairpin.fa'
    mature_url = 'https://www.mirbase.org/download/mature.fa'

    response = requests.get(hairpin_url)
    if response.status_code == 200:
        with open(f'hairpin.fa', 'wb') as f:
            f.write(response.content)
    else:
        print(f'Error when retrieving hairpin.fa\nError: {response.status_code}')
        
    response = requests.get(mature_url)
    if response.status_code == 200:
        with open(f'mature.fa', 'wb') as f:
            f.write(response.content)
    else:
        print(f'Error when retrieving mature.fa\nError: {response.status_code}')

# TODO figure out how to filter out non-GATCU chracters. 
def filter_miRNA_refs():
    # filter for mature human and replace non ACGTU letter with N
    
    with open("mature_hsa.fa", "w") as hsa_file, open("mature_other.fa", "w") as other_file:
        for record in SeqIO.parse("mature.fa", "fasta"):
            temp = record.name
            if "hsa" in temp:
                #sequence = re.sub('[^GATCU]', "N", str(record.seq).upper())
                #record.seq = sequence
                SeqIO.write(record, hsa_file, "fasta")
            else:
                #sequence = re.sub('[^GATCU]', "N", str(record.seq).upper())
                #record.seq = sequence
                SeqIO.write(record, other_file, "fasta")
        

    with open("hairpin_hsa.fa", "w") as hsa_file:
        for record in SeqIO.parse("hairpin.fa", "fasta"):
            temp = record.name
            if "hsa" in temp:
                #sequence = re.sub('[^GATCU]', "N", str(record.seq).upper())
                #record.seq = sequence
                SeqIO.write(record, hsa_file, "fasta")

def main():
    print("Generate the miRNA reference files, separating human (hsa) from other organisms.")
    print("Files will be generated in the directory you are running this script.")

    download_miRBase()

    filter_miRNA_refs()

    print("Completed")

if __name__ == "__main__":
    main()