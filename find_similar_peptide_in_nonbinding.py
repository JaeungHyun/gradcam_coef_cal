import pandas as pd
import re
from tqdm import tqdm
import sys
from util import call_group_list, find_group


allele = sys.argv[1]

df = pd.read_pickle('/home/jaeung/970evo/MHC/DeepNeo_new_testset.pkl')
del df['matrix']
df['length'] = df['Peptide seq'].map(lambda x: len(x))
df = df[df['length'] == 9]

df3 = df[df['answer'] == 0]
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
    try:
        binder = pd.read_csv(f'{allele} {target_list[j]} p2,9.csv')

        df2 = df[df['allele'].isin(target)]
        df4 = df3[df3['allele'].isin(target)]
        #pep_list = df2['Peptide seq'].unique().tolist()
        pep_list = binder['Peptide seq'].unique().tolist()
        pep_list.sort()
        print(pep_list
              )
        total_df = []
        for pep1 in pep_list:
            globals()[f'{pep1}_similar'] = []
            tmp = ''
            for i, p in enumerate(pep1):
                # if i == check_position[j]:
                if i == 1 or i == 8:
                    tmp += '[A-Z]'
                else:
                    tmp += p
            p = re.compile(tmp)

            for pep2 in df4['Peptide seq'].unique():
                if p.search(pep2):
                    globals()[f'{pep1}_similar'].append(pep2)

        target = []
        for pep1 in pep_list:
            #if len(globals()[f'{pep1}_similar']) > 1:
            target.append(pep1)

        target_df = []
        for t in target:
            tmp = df4[df4['Peptide seq'].isin(globals()[f'{t}_similar'])]
            target_df.append(tmp)
    
        result_df = []
        for t_df in target_df:
            t_df = t_df.sort_values('Peptide seq')
            t_df.reset_index(drop=True, inplace=True)
            t_df['difference_position'] = check_position[j]+1
            t_df['Hydrophobicity'] = ''
            t_df['Bulkiness'] = ''
            t_df['Polarity'] = ''
            t_df['Charge'] = ''
            t_df['MW'] = ''

            for i in range(len(t_df)):
                if t_df['difference_position'][i]:
                    t_df.at[i, 'Hydrophobicity'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check_position[j]], 'Hydrophobicity'].values[0]
                    t_df.at[i, 'Bulkiness'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check_position[j]], 'Bulkiness'].values[0]
                    t_df.at[i, 'Charge'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check_position[j]], 'Charge'].values[0]
                    t_df.at[i, 'MW'] = aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check_position[j]], 'MW'].values[0]
                    t_df.at[i, 'Polarity'] = \
                        aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][check_position[j]], 'Polarity'].values[0]
            result_df.append(t_df)

        for _, t_df in enumerate(result_df):
            t_df['group'] = t_df['allele'].apply(find_group)

        for tdf in result_df:
            total_df.append(tdf)
        try:
            pd.concat(total_df).drop_duplicates().to_csv(f'{allele} group{j + 1} nonbinding p2,9.csv')
        except:
            pass
    except:
        print(allele, target_list[j])