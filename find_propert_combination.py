import pandas as pd
import pickle
import numpy as np
import gc
from util import *
import ray
from itertools import product


def check_combi(hla, pep, mode):
    target1 = aa_property[aa_property[mode] == 1]['aa'].tolist()
    target2 = aa_property[aa_property[mode] == 0]['aa'].tolist()

    if hla in target1 and pep in target1:
        return 0
    elif hla in target1 and pep in target2:
        return 1
    elif hla in target2 and pep in target1:
        return 2
    elif hla in target2 and pep in target2:
        return 3


@ray.remote
def find_property(df, target_group, binder, hla, allele, target, mode):
    result = {}
    result[allele] = []
    df = df[df['allele'].isin(target_group)]
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        new_array = np.zeros((hla_len, 9))
        for h, i in enumerate(hla[allele]):
            for w, j in enumerate(pepseq):
                if check_combi(i, j, mode) == target:
                    new_array[h][w] = binder[allele][num][h][w]

        result[allele].append(new_array)

    return result


def load_gradcam_result():
    with open('/home/jaeung/Research/MHC/ms+ba_short_hla_gradcam_result.pkl', 'rb') as f:
        return pickle.load(f)


def load_short_hla():
    hla_b_prot = pd.read_csv('/home/jaeung/Research/MHC/HLA_B_prot.txt', sep='\t', header = None)
    hla_a_prot = pd.read_csv('/home/jaeung/Research/MHC/HLA_A_prot.txt', sep='\t', header = None)
    hla_c_prot = pd.read_csv('/home/jaeung/Research/MHC/HLA_C_prot.txt', sep='\t', header = None)

    hla_a_prot[1] = hla_a_prot[1].map(lambda x: x[24:-65])
    hla_c_prot[1] = hla_c_prot[1].map(lambda x: x[4:-66])
    hla_b_prot[1] = hla_b_prot[1].map(lambda x: x[12:-62])
    hla_prot = pd.concat([hla_a_prot, hla_b_prot, hla_c_prot], axis = 0)

    hla = {}
    for line in hla_prot.to_numpy():
        hla[line[0]] = line[1]
    return hla


def load_pep_seq():
    with open('/home/jaeung/Research/MHC/Short_HLA_seq_training_data.pkl', 'rb') as f:
        df = pickle.load(f)

    del df['matrix'], df['sequence']
    gc.collect()

    df['length'] = df['Peptide seq'].map(lambda x: len(x))
    df = df[df['length'] == 9]
    return df


if __name__ == "main":
    ray.init(address='auto')

    p9_binder, _, _, _ = load_gradcam_result()
    df = load_pep_seq()
    hla = load_short_hla()

    aa_property = pd.read_excel('Amino_acid_property.xlsx')
    aa_property['hydro'] = aa_property['Hydrophobicity'].map(lambda x: 1 if x >= 0 else 0)
    aa_property['bulky'] = aa_property['Bulkiness'].map(lambda x: 1 if x >= 15.4 else 0)  # 평균값이 15.367500000000001

    p9_binder_id = ray.put(p9_binder)

    hla_id = ray.put(hla)
    df_id = ray.put(df)

    item = [['HLA-A', 'HLA-B', 'HLA-C'], ['bulky', 'hydro'], [0, 1, 2, 3]]

    for allele, mode, target in list(product(*item)):
        if allele == 'HLA-A':
            hla_len = 276
        elif allele == 'HLA-B':
            hla_len = 252
        else:
            hla_len = 268

        target_list, group_list = call_group_list(allele)
        total_g = []
        for g in group_list:
            total_g.extend(g)

        print('Start Ray')
        result = ray.get([find_property.remote(df_id, total_g, p9_binder_id, hla_id, allele, target, mode)
                          for allele in total_g])

        print('Saving Result')
        with open(f'/home/jaeung/Research/MHC/{allele}_{mode}_gradcam_result.pkl', 'wb') as f:
            pickle.dump(result, f)

