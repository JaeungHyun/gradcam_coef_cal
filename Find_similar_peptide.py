import pandas as pd
import pickle
import numpy as np
import re
from tqdm import tqdm
import sys
from util import call_group_list


def find_group(allele):
    if allele in group1_a:
        return 'group1'
    elif allele in group1_a:
        return 'group1'
    elif allele in group2_a:
        return 'group2'
    elif allele in group3_a:
        return 'group3'
    elif allele in group4_a:
        return 'group4'
    elif allele in group5_a:
        return 'group5'
    elif allele in group6_a:
        return 'group6'
    elif allele in group7_a:
        return 'group7'


    elif allele in group1_b:
        return 'group1'
    elif allele in group2_b:
        return 'group2'
    elif allele in group3_b:
        return 'group3'
    elif allele in group4_b:
        return 'group4'
    elif allele in group5_b:
        return 'group5'
    elif allele in group6_b:
        return 'group6'
    elif allele in group7_b:
        return 'group7'
    elif allele in group7_b:
        return 'group8'
    elif allele in group9_b:
        return 'group9'
    elif allele in group10_b:
        return 'group10'
    elif allele in group11_b:
        return 'group11'
    elif allele in group12_b:
        return 'group12'
    elif allele in group13_b:
        return 'group13'
    elif allele in group14_b:
        return 'group14'
    elif allele in group15_b:
        return 'group15'
    elif allele in group16_b:
        return 'group16'
    elif allele in group17_b:
        return 'group17'

    elif allele in group1_c:
        return 'group1'
    elif allele in group2_c:
        return 'group2'
    elif allele in group3_c:
        return 'group3'
    elif allele in group4_c:
        return 'group4'
    elif allele in group5_c:
        return 'group5'


allele = sys.argv[1]

df = pd.read_pickle('/home/jaeung/Research/MHC/Model data/2021_testset+before_dataset.pkl')
del df['matrix']
df['length'] = df['Peptide seq'].map(lambda x: len(x))
df = df[df['length'] == 9]
df = df[df['answer'] == 1]

group1_a = ['HLA-A-2403', 'HLA-A-2402', 'HLA-A-2413', 'HLA-A-2301', 'HLA-A-2406', 'HLA-A-2407']
group2_a = ['HLA-A-3303', 'HLA-A-3301', 'HLA-A-6801', 'HLA-A-6601', 'HLA-A-3401', 'HLA-A-6602',
            'HLA-A-3101', 'HLA-A-7401']
group3_a = ['HLA-A-3001', 'HLA-A-0301', 'HLA-A-1101', 'HLA-A-1102', 'HLA-A-6812']
group4_a = ['HLA-A-6802', 'HLA-A-6901']
group5_a = ['HLA-A-0205', 'HLA-A-0206', 'HLA-A-0217', 'HLA-A-0216', 'HLA-A-0212', 'HLA-A-0219',
            'HLA-A-0207', 'HLA-A-0203', 'HLA-A-0201', 'HLA-A-0211', 'HLA-A-0204', 'HLA-A-0202']
group6_a = ['HLA-A-2601', 'HLA-A-2501', 'HLA-A-2608', 'HLA-A-2603', 'HLA-A-2602']
group7_a = ['HLA-A-0103', 'HLA-A-0101', 'HLA-A-2902', 'HLA-A-3002', 'HLA-A-3601', 'HLA-A-8001']
group1_b = ['HLA-B-5301', 'HLA-B-3501', 'HLA-B-3507', 'HLA-B-3508', 'HLA-B-1511']
group2_b = ['HLA-B-0704', 'HLA-B-0702', 'HLA-B-4201', 'HLA-B-3502', 'HLA-B-3503', 'HLA-B-3504', 'HLA-B-3506', ]
group3_b = ['HLA-B-8101', 'HLA-B-4202', ]
group4_b = ['HLA-B-5401', 'HLA-B-5501', ]
group5_b = ['HLA-B-5502', 'HLA-B-5601']
group6_b = ['HLA-B-5101', 'HLA-B-5108', 'HLA-B-7301', 'HLA-B-3906', ]
group7_b = ['HLA-B-2710', 'HLA-B-2702', 'HLA-B-2701', 'HLA-B-2704', 'HLA-B-2703', 'HLA-B-2705', 'HLA-B-2708',
            'HLA-B-2707', 'HLA-B-2706', ]
group8_b = ['HLA-B-3905', 'HLA-B-3901', 'HLA-B-3801', 'HLA-B-3802', 'HLA-B-1509', 'HLA-B-1510', ]
group9_b = ['HLA-B-3924', 'HLA-B-1402', 'HLA-B-1403', ]
group10_b = ['HLA-B-2709', 'HLA-B-3909', ]
group11_b = ['HLA-B-4901', 'HLA-B-5001', 'HLA-B-4006', 'HLA-B-4101', 'HLA-B-4501', ]
group12_b = ['HLA-B-1803', 'HLA-B-1801', 'HLA-B-4402', 'HLA-B-4403', 'HLA-B-4427', 'HLA-B-4428', ]
group13_b = ['HLA-B-4102', 'HLA-B-4104', 'HLA-B-4103', 'HLA-B-4409', 'HLA-B-4002', 'HLA-B-4001', ]
group14_b = ['HLA-B-1508', 'HLA-B-1501', 'HLA-B-1503', 'HLA-B-1502', 'HLA-B-4601', ]
group15_b = ['HLA-B-5703', 'HLA-B-5701', 'HLA-B-5801', 'HLA-B-5802', 'HLA-B-1517', ]
group16_b = ['HLA-B-5201', 'HLA-B-1302', ]
group17_b = ['HLA-B-0803', 'HLA-B-0802', ]
group1_c = ['HLA-C-0401', 'HLA-C-0501', 'HLA-C-0403', 'HLA-C-0802', ]
group2_c = ['HLA-C-1402', 'HLA-C-1402', ]
group3_c = ['HLA-C-0704', 'HLA-C-0702', 'HLA-C-0602', 'HLA-C-0701', ]
group4_c = ['HLA-C-1502', 'HLA-C-1505', ]
group5_c = ['HLA-C-1701', 'HLA-C-0801', 'HLA-C-0304', 'HLA-C-0303', 'HLA-C-1202',
            'HLA-C-0202', 'HLA-C-1203', 'HLA-C-1601', 'HLA-C-0302', ]

if allele == 'HLA-A':
    check_position = [8] * 7
elif allele == 'HLA-B':
    check_position = [1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1]
else:
    check_position = [1, 1, 1, 1, 1]

target_list, group_list = call_group_list(allele)

for j, target in enumerate(group_list):
    df2 = df[df['allele'].isin(target)]
    pep_list = df2['Peptide seq'].unique().tolist()
    pep_list.sort()
    tmp = ''
    for i, pep1 in tqdm(enumerate(pep_list)):
        globals()[f'{pep1}_similar'] = []
        tmp = ''
        for i, p in enumerate(pep1):
            # if i == 1 or i == 8:
            if i == check_position[j]:
                tmp += '[A-Z]'
            else:
                tmp += p
        p = re.compile(tmp)

        for pep2 in pep_list:
            if p.search(pep2):
                globals()[f'{pep1}_similar'].append(pep2)

    count = []
    target = []
    for i, pep1 in tqdm(enumerate(pep_list)):
        if len(globals()[f'{pep1}_similar']) > 1:
            count.append(i)
            target.append(pep1)
        else:
            pass

    target_df = []
    for i, t in tqdm(enumerate(target)):
        tmp = df2[df2['Peptide seq'].isin(globals()[f'{t}_similar'])]
        target_df.append(tmp)
    aa_property = pd.read_excel('Amino_acid_property.xlsx')

    result_df = []
    for t_df in target_df:
        t_df = t_df.sort_values('Peptide seq')
        t_df.reset_index(drop=True, inplace=True)
        t_df['difference_position_2'] = ''
        t_df['difference_position_9'] = ''
        t_df['Hydrophobicity_2'] = ''
        t_df['Bulkiness_2'] = ''
        t_df['Polarity_2'] = ''
        t_df['Charge_2'] = ''
        t_df['MW_2'] = ''
        t_df['Hydrophobicity_9'] = ''
        t_df['Bulkiness_9'] = ''
        t_df['Polarity_9'] = ''
        t_df['Charge_9'] = ''
        t_df['MW_9'] = ''

        for i in range(1, len(t_df)):
            if t_df['Peptide seq'][0][1] != t_df['Peptide seq'][i][1]:
                t_df.at[i, 'difference_position_2'] = 2
            if t_df['Peptide seq'][0][8] != t_df['Peptide seq'][i][8]:
                t_df.at[i, 'difference_position_9'] = 9
            # 기준점 말고 다른 부분에도 데이터 추가
        for i in range(len(t_df)):
            if t_df.at[i, 'difference_position_2'] == '':
                t_df.at[i, 'difference_position_2'] = t_df['difference_position_2'].unique()[-1]
            if t_df.at[i, 'difference_position_9'] == '':
                t_df.at[i, 'difference_position_9'] = t_df['difference_position_9'].unique()[-1]
        for i in range(1, len(t_df)):
            if t_df.at[i, 'Peptide seq'] == t_df.at[0, 'Peptide seq']:
                t_df.at[i, 'difference_position_2'] = t_df.at[0, 'difference_position_2']
                t_df.at[i, 'difference_position_9'] = t_df.at[0, 'difference_position_9']
        for i in range(len(t_df)):
            if t_df['difference_position_2'][i]:
                t_df.at[i, 'Hydrophobicity_2'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][1], 'Hydrophobicity'].values[0]
                t_df.at[i, 'Bulkiness_2'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][1], 'Bulkiness'].values[0]
                t_df.at[i, 'Charge_2'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][1], 'Charge'].values[0]
                t_df.at[i, 'MW_2'] = aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][1], 'MW'].values[0]
                t_df.at[i, 'Polarity_2'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][1], 'Polarity'].values[0]
            if t_df['difference_position_9'][i]:
                t_df.at[i, 'Hydrophobicity_9'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][8], 'Hydrophobicity'].values[0]
                t_df.at[i, 'Bulkiness_9'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][8], 'Bulkiness'].values[0]
                t_df.at[i, 'Charge_9'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][8], 'Charge'].values[0]
                t_df.at[i, 'MW_9'] = aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][8], 'MW'].values[0]
                t_df.at[i, 'Polarity_9'] = \
                    aa_property.loc[aa_property['aa'] == t_df['Peptide seq'][i][8], 'Polarity'].values[0]
        result_df.append(t_df)

    for t_df in result_df:
        t_df['difference_position'] = ''
        t_df['Hydrophobicity'] = ''
        t_df['Bulkiness'] = ''
        t_df['Polarity'] = ''
        t_df['Charge'] = ''
        t_df['MW'] = ''

        for i in range(len(t_df)):
            if t_df['difference_position_2'][i] != '' and t_df['difference_position_9'][i] != '':
                t_df.at[i, 'difference_position'] = str(t_df['difference_position_2'][i]) + ';' + str(
                    t_df['difference_position_9'][i])
                t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_2'][i]) + ';' + str(t_df['Hydrophobicity_9'][i])
                t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_2'][i]) + ';' + str(t_df['Bulkiness_9'][i])
                t_df.at[i, 'Polarity'] = str(t_df['Polarity_2'][i]) + ';' + str(t_df['Polarity_9'][i])
                t_df.at[i, 'Charge'] = str(t_df['Charge_2'][i]) + ';' + str(t_df['Charge_9'][i])
                t_df.at[i, 'MW'] = str(t_df['MW_2'][i]) + ';' + str(t_df['MW_9'][i])

            elif t_df['difference_position_2'][i] != '' and t_df['difference_position_9'][i] == '':
                t_df.at[i, 'difference_position'] = str(t_df['difference_position_2'][i])
                t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_2'][i])
                t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_2'][i])
                t_df.at[i, 'Polarity'] = str(t_df['Polarity_2'][i])
                t_df.at[i, 'Charge'] = str(t_df['Charge_2'][i])
                t_df.at[i, 'MW'] = str(t_df['MW_2'][i])

            elif t_df['difference_position_2'][i] == '' and t_df['difference_position_9'][i] != '':
                t_df.at[i, 'difference_position'] = str(t_df['difference_position_9'][i])
                t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_9'][i])
                t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_9'][i])
                t_df.at[i, 'Polarity'] = str(t_df['Polarity_9'][i])
                t_df.at[i, 'Charge'] = str(t_df['Charge_9'][i])
                t_df.at[i, 'MW'] = str(t_df['MW_9'][i])
        for i in range(len(t_df)):
            if t_df['difference_position_2'][i] == '' and t_df['difference_position_9'][i] == '':
                if t_df['difference_position_2'].unique()[-1] != '' and t_df['difference_position_9'].unique()[-1] != '':
                    t_df.at[i, 'difference_position'] = str(t_df['difference_position_2'][i]) + ';' + str(
                        t_df['difference_position_9'][i])
                    t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_2'][i]) + ';' + str(
                        t_df['Hydrophobicity_9'][i])
                    t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_2'][i]) + ';' + str(t_df['Bulkiness_9'][i])
                    t_df.at[i, 'Polarity'] = str(t_df['Polarity_2'][i]) + ';' + str(t_df['Polarity_9'][i])
                    t_df.at[i, 'Charge'] = str(t_df['Charge_2'][i]) + ';' + str(t_df['Charge_9'][i])
                    t_df.at[i, 'MW'] = str(t_df['MW_2'][i]) + ';' + str(t_df['MW_9'][i])
                elif t_df['difference_position_2'].unique()[-1] != '' and t_df['difference_position_9'].unique()[-1] == '':
                    t_df.at[i, 'difference_position'] = str(t_df['difference_position_2'][i])
                    t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_2'][i])
                    t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_2'][i])
                    t_df.at[i, 'Polarity'] = str(t_df['Polarity_2'][i])
                    t_df.at[i, 'Charge'] = str(t_df['Charge_2'][i])
                    t_df.at[i, 'MW'] = str(t_df['MW_2'][i])
                elif t_df['difference_position_2'].unique()[-1] == '' and t_df['difference_position_9'].unique()[-1] != '':
                    t_df.at[i, 'difference_position'] = str(t_df['difference_position_9'][i])
                    t_df.at[i, 'Hydrophobicity'] = str(t_df['Hydrophobicity_9'][i])
                    t_df.at[i, 'Bulkiness'] = str(t_df['Bulkiness_9'][i])
                    t_df.at[i, 'Polarity'] = str(t_df['Polarity_9'][i])
                    t_df.at[i, 'Charge'] = str(t_df['Charge_9'][i])
                    t_df.at[i, 'MW'] = str(t_df['MW_9'][i])

    for i, t_df in enumerate(result_df):
        result_df[i] = t_df[['allele',
                             'answer',
                             'Peptide seq',
                             'length',
                             'difference_position',
                             'Hydrophobicity',
                             'Bulkiness',
                             'Polarity',
                             'Charge',
                             'MW']]

    for i, t_df in enumerate(result_df):
        t_df['group'] = t_df['allele'].apply(find_group)
    tt = []
    for tdf in result_df:
        tt.append(tdf)
    try:
        pd.concat(tt).drop_duplicates().to_csv(f'{allele} group{j + 1} except p2,p9.csv')
    except:
        print(target)
