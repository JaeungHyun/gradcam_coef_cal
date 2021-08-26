import pandas as pd
import re
from tqdm import tqdm
import sys
from util import call_group_list, find_group


allele = sys.argv[1]

df = pd.read_pickle('/home/jaeung/Research/MHC/Model data/2021_testset+before_dataset.pkl')
del df['matrix']
df['length'] = df['Peptide seq'].map(lambda x: len(x))
df = df[df['length'] == 9]
df = df[df['answer'] == 1]
aa_property = pd.read_excel('Amino_acid_property.xlsx')


if allele == 'HLA-A':
    check_position = [8] * 7
elif allele == 'HLA-B':
    check_position = [1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1]
else:
    check_position = [1, 1, 1, 1, 1]

target_list, group_list = call_group_list(allele)

for j, target in enumerate(tqdm(group_list)):
    df2 = df[df['allele'].isin(target)]
    pep_list = df2['Peptide seq'].unique().tolist()
    pep_list.sort()
    check_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    check_list.pop(check_position[j])
    total_df = []
    for check in check_list:
        for pep1 in pep_list:
            globals()[f'{pep1}_similar'] = []
            tmp = ''
            for i, p in enumerate(pep1):
                if i == check:
                    tmp += '[A-Z]'
                else:
                    tmp += p
            p = re.compile(tmp)

            for pep2 in pep_list:
                if p.search(pep2):
                    globals()[f'{pep1}_similar'].append(pep2)

        target = []
        for pep1 in pep_list:
            if len(globals()[f'{pep1}_similar']) > 1:
                target.append(pep1)
            else:
                pass

        target_df = []
        for t in target:
            tmp = df2[df2['Peptide seq'].isin(globals()[f'{t}_similar'])]
            target_df.append(tmp)

        result_df = []
        for t_df in target_df:
            t_df = t_df.sort_values('Peptide seq')
            t_df.reset_index(drop=True, inplace=True)
            t_df['difference_position'] = check+1
            t_df['Hydrophobicity'] = ''
            t_df['Bulkiness'] = ''
            t_df['Polarity'] = ''
            t_df['Charge'] = ''
            t_df['MW'] = ''

            for i in range(len(t_df)):
                if t_df['difference_position'][i]:
                    t_df.at[i, 'Hydrophobicity'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check], 'Hydrophobicity'].values[0]
                    t_df.at[i, 'Bulkiness'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check], 'Bulkiness'].values[0]
                    t_df.at[i, 'Charge'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check], 'Charge'].values[0]
                    t_df.at[i, 'MW'] = aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check], 'MW'].values[0]
                    t_df.at[i, 'Polarity'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check], 'Polarity'].values[0]
            result_df.append(t_df)

        for _, t_df in enumerate(result_df):
            t_df['group'] = t_df['allele'].apply(find_group)

        for tdf in result_df:
            total_df.append(tdf)
    try:
        pd.concat(total_df).drop_duplicates().to_csv(f'{allele} group{j + 1} Position total.csv')
    except:
        pass
