import numpy as np
from util import *
import ray
from itertools import product
import sys


def check_combi(pep, mode):
    target1 = aa_property[aa_property[mode] == 1]['aa'].tolist()
    target2 = aa_property[aa_property[mode] == 0]['aa'].tolist()

    if pep in target1:
        return 0
    elif pep in target2:
        return 1


@ray.remote
def find_property(df, target_group, binder, allele, target, mode, p):
    result = {}
    result[allele] = []
    df = df[df['allele'].isin(target_group)]  # HLA-A,B,C 각각 가져오는 부분
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        if check_combi(pepseq[p], mode) == target:
            result[allele].append(binder[allele][num][:, p])
            # new_array[:, w] = binder[allele][num][:, w]

        # result[allele].append(new_array)

    return result


def find_property_value(aa, mode):
    if mode == 'polar':
        v = 4
    elif mode == 'bulky':
        v = 3
    elif mode == 'hydro':
        v = 2
    elif mode == 'Charge':
        v = 5
    elif mode == 'MW':
        v = 6
    # print(v)
    value = aa_property.loc[aa_property['aa'] == aa.upper()].values[0][v]
    # print(value)
    try:
        return value
    except:
        print(aa.upper())


@ray.remote
def find_property2(df, target_group, binder, allele, mode, p):
    cp_value = {}
    cor_result = {}

    cor_result[allele] = []
    cp_value[allele] = []
    df = df[df['allele'].isin(target_group)]  # HLA-A,B,C 각각 가져오는 부분
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq'].to_numpy()):
        # print(p)
        # print(num)
        cor_result[allele].append(binder[allele][num][:, p])
        cp_value[allele].append(find_property_value(pepseq[p], mode))

    return cp_value, cor_result


@ray.remote
def find_property3(df, target_group, binder, allele, mode, p):
    cp_value = {}
    cp_value[allele] = []
    df = df[df['allele'].isin(target_group)]  # HLA-A,B,C 각각 가져오는 부분
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        #print(pepseq[p])
        cp_value[allele].append(find_property_value(pepseq[p], mode))
    return cp_value


# if __name__ == "main":


try:
    ray.init(address='auto')
except:
    ray.init(dashboard_host='0.0.0.0',
             log_to_driver=False,
             _plasma_directory="/home/jaeung/tmp"
             )

p9_binder = load_gradcam_result(sys.argv[4])
p9_binder_id = ray.put(p9_binder)
df = load_pep_seq()
hla = load_short_hla()
hla_id = ray.put(hla)
df_id = ray.put(df)
del p9_binder, df, hla

aa_property = pd.read_excel('Amino_acid_property.xlsx')

item = [[sys.argv[1]], ['polar','hydro','bulky', 'MW', 'Charge']]

for allele, mode in list(product(*item)):
    print(allele, mode)
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
    for p in range(9):
        if sys.argv[2] == 'random':
            result = ray.get([find_property3.remote(df_id, total_g, p9_binder_id, allele, mode, p)
                              for allele in total_g])
            print('Saving Result')
            with open(
                    f'/home/jaeung/960evo/result/{allele}_random_protein_{mode}_position_{p + 1}_gradcam_result_with_cp_value.pkl',
                    'wb') as f:
                pickle.dump(result, f)
        elif sys.argv[4] == 'natural':
            result = ray.get([find_property3.remote(df_id, total_g, p9_binder_id, allele, mode, p)
                              for allele in total_g])
            print('Saving Result')
            with open(
                    f'/home/jaeung/960evo/result/{allele}_natural_protein_{mode}_position_{p + 1}_gradcam_result_with_cp_value.pkl',
                    'wb') as f:
                pickle.dump(result, f)


