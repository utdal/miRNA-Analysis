import pandas as pd
import requests

class GetmiRNATargets:
    def __init__(self, miRNA_list):
        self.miRNA_list = miRNA_list

    # Function to get the data about the miRNA from ENCORI
    # https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all&cellType=all'

    def retrive_miRNA_csvs(self,miRNA):
        url = f'https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA={miRNA}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all&cellType=all'
        response = requests.get(url)
        if response.status_code == 200:
            with open(f'{miRNA}_targets.csv', 'wb') as f:
                f.write(response.content)
        else:
            print(f'Error: {response.status_code}')

    def create_df_of_all_data(self):

        miRNA_targets = pd.DataFrame()

        # Get the csvs from the ENCORI database
        # TODO consider adding options for restrictions on info pulled.
        for index,row in self.miRNA_list.iterrows():
            self.retrive_miRNA_csvs(row['miRNA'])

            # write the data from the csv into miRNA_targets DataFrame with all miRNA and their targets

