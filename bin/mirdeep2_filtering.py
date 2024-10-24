import pandas as pd

# Get novel miRNAs detected data
m1novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM1_S5.tsv', sep='\t', header=0)
m2novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM2.tsv', sep='\t', header=0)
m4novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM4.tsv', sep='\t', header=0)
f1novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelf1.tsv', sep='\t', header=0)
f2novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelf2.tsv', sep='\t', header=0)
f3novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelf3.tsv', sep='\t', header=0)
f4novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelf4.tsv', sep='\t', header=0)
m20novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM20.csv', sep='\t', header=0)
m10novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM10.csv', sep='\t', header=0)
m21novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM21.csv', sep='\t', header=0)
m22novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM22.csv', sep='\t', header=0)
m23novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM23.csv', sep='\t', header=0)
m24novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM24.csv', sep='\t', header=0)
m26novel = pd.read_csv('{base_dir}/Analysis_Scripts/novelM26.csv', sep='\t', header=0)
x80T2Rnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel80T2R.csv', sep='\t', header=0)
x83T8Rnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel83T8R.csv', sep='\t', header=0)
x84T4Lnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel84T4L.csv', sep='\t', header=0)
x84T5Lnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel84T5L.csv', sep='\t', header=0)
x85T3Lnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel85T3L.csv', sep='\t', header=0)
x85T3Rnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel85T3R.csv', sep='\t', header=0)
x85T5Rnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel85T5R.csv', sep='\t', header=0)
x91T5Rnovel = pd.read_csv('{base_dir}/Analysis_Scripts/novel91T5R.csv', sep='\t', header=0)

# add a column for the sample name
m1novel.insert(0,'Sample Name', 'M1_S5')
m2novel.insert(0,'Sample Name', 'M2_S6')
m4novel.insert(0,'Sample Name', 'M4_S7')
f1novel.insert(0,'Sample Name', 'F1_S1')
f2novel.insert(0,'Sample Name', 'F2_S2')
f3novel.insert(0,'Sample Name', 'F3_S3')
f4novel.insert(0,'Sample Name', 'F4_S4')
m20novel.insert(0,'Sample Name', 'M20_S9')
m10novel.insert(0,'Sample Name', 'M10_S8')
m21novel.insert(0,'Sample Name', 'M21_S10')
m22novel.insert(0,'Sample Name', 'M22_S11')
m23novel.insert(0,'Sample Name', 'M23_S12')
m24novel.insert(0,'Sample Name', 'M24_S13')
m26novel.insert(0,'Sample Name', 'M26_S14')
x80T2Rnovel.insert(0,'Sample Name', '80T2R_S17')
x83T8Rnovel.insert(0,'Sample Name', '83T8R_S18')
x84T4Lnovel.insert(0,'Sample Name', '84T4L_S15')
x84T5Lnovel.insert(0,'Sample Name', '84T5L_S19')
x85T3Lnovel.insert(0,'Sample Name', '85T3L_S21')
x85T3Rnovel.insert(0,'Sample Name', '85T3R_S16')
x85T5Rnovel.insert(0,'Sample Name', '85T5R_S22')
x91T5Rnovel.insert(0,'Sample Name', '91T5R_S20')

# Remove rows where the miRDeep2 score is less than 4. 
# 4 is supposedly the best liked cutoff score
m1novel = m1novel.drop(m1novel[m1novel['miRDeep2 score'] < 4].index)
m2novel = m2novel.drop(m2novel[m2novel['miRDeep2 score'] < 4].index)
m4novel = m4novel.drop(m4novel[m4novel['miRDeep2 score'] < 4].index)
f1novel = f1novel.drop(f1novel[f1novel['miRDeep2 score'] < 4].index)
f2novel = f2novel.drop(f2novel[f2novel['miRDeep2 score'] < 4].index)
f3novel = f3novel.drop(f3novel[f3novel['miRDeep2 score'] < 4].index)
f4novel = f4novel.drop(f4novel[f4novel['miRDeep2 score'] < 4].index)
m20novel = m20novel.drop(m20novel[m20novel['miRDeep2 score'] < 4].index)
m10novel = m10novel.drop(m10novel[m10novel['miRDeep2 score'] < 4].index)
m21novel = m21novel.drop(m21novel[m21novel['miRDeep2 score'] < 4].index)
m22novel = m22novel.drop(m22novel[m22novel['miRDeep2 score'] < 4].index)
m23novel = m23novel.drop(m23novel[m23novel['miRDeep2 score'] < 4].index)
m24novel = m24novel.drop(m24novel[m24novel['miRDeep2 score'] < 4].index)
m26novel = m26novel.drop(m26novel[m26novel['miRDeep2 score'] < 4].index)
x80T2Rnovel = x80T2Rnovel.drop(x80T2Rnovel[x80T2Rnovel['miRDeep2 score'] < 4].index)
x83T8Rnovel = x83T8Rnovel.drop(x83T8Rnovel[x83T8Rnovel['miRDeep2 score'] < 4].index)
x84T4Lnovel = x84T4Lnovel.drop(x84T4Lnovel[x84T4Lnovel['miRDeep2 score'] < 4].index)
x84T5Lnovel = x84T5Lnovel.drop(x84T5Lnovel[x84T5Lnovel['miRDeep2 score'] < 4].index)
x85T3Lnovel = x85T3Lnovel.drop(x85T3Lnovel[x85T3Lnovel['miRDeep2 score'] < 4].index)
x85T3Rnovel = x85T3Rnovel.drop(x85T3Rnovel[x85T3Rnovel['miRDeep2 score'] < 4].index)
x85T5Rnovel = x85T5Rnovel.drop(x85T5Rnovel[x85T5Rnovel['miRDeep2 score'] < 4].index)
x91T5Rnovel = x91T5Rnovel.drop(x91T5Rnovel[x91T5Rnovel['miRDeep2 score'] < 4].index)

# remove columns I don't want
# remove rfam alert, mirbase mirna, ucsc browser, ncbi blast
m1novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m2novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m4novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
f1novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
f2novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
f3novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
f4novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m20novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m10novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m21novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m22novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m23novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m24novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
m26novel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x80T2Rnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x83T8Rnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x84T4Lnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x84T5Lnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x85T3Lnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x85T3Rnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x85T5Rnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)
x91T5Rnovel.drop(['rfam alert','miRBase miRNA','UCSC browser','NCBI blastn'], axis=1)

# Make one dataframe with all samples. 
# There is still column that defines which sample it is from
all_samples = pd.concat([m1novel,m2novel,m4novel,f1novel,f2novel,f3novel,f4novel,m20novel,m10novel,m21novel,m22novel,m23novel,m24novel,m26novel,x80T2Rnovel,x83T8Rnovel,x84T4Lnovel,x84T5Lnovel,x85T3Lnovel,x85T3Rnovel,x85T5Rnovel,x91T5Rnovel])

class FilterMiRDeep2Results:
    def __init__(self, mirdeep_res_loc):
        self.mirdeep_res_loc = mirdeep_res_loc

    # Function for if the precursor sequences are amost the same?
    def sim_hairpin(self, hairpin1, hairpin2):
        if hairpin1 in hairpin2 or hairpin2 in hairpin1:
            return True
        else:
            return False
        
    # Create a dataframe with the similar precursors in order
    columns = ['Sample Name', 'provisional id', 'miRDeep2 score',
        'estimated probability that the miRNA candidate is a true positive',
        'total read count', 'mature read count',
        'loop read count', 'star read count', 'significant randfold p-value',
        'miRBase miRNA', 'example miRBase miRNA with the same seed',
        'consensus mature sequence',
        'consensus star sequence', 'consensus precursor sequence',
        'precursor coordinate']

    similar_hairpin = pd.DataFrame(columns=columns)
    only_one = pd.DataFrame(columns=columns)

    similar_hairpin_rows = []
    only_one_rows = []

    visited_indices = set()

    for test_index,test_row in all_samples.iterrows():
        if test_index in visited_indices:
            continue

        has_sim = False
        for index,row in all_samples.iterrows():
            if index != test_index and index not in visited_indices and sim_hairpin(test_row['consensus precursor sequence'], row['consensus precursor sequence']):
                # copy row into similar_hairpin + delete row
                has_sim = True
                similar_hairpin.loc[len(similar_hairpin)] = row
                visited_indices.add(index)

        if has_sim is True:
            # write test_row to the similarl_hairpin
            similar_hairpin.loc[len(similar_hairpin)] = test_row
        else:
            # write test_row to only_one
            only_one.loc[len(only_one)] = test_row

    similar_hairpin.to_csv('{base_dir}/Analysis_Scripts/similar_hairpin_all_samples.csv', sep='\t')

    only_one.to_csv('{base_dir}/Analysis_Scripts/no_repeat_across_all_samples.csv', sep='\t')