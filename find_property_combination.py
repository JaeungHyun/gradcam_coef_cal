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
    df = df[df['allele'].isin(target_group)] # HLA-A,B,C 각각 가져오는 부분
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        #new_array = np.zeros((hla_len, 9))
        #for h, i in enumerate(hla[allele]):
        #for w, j in enumerate(pepseq):
        if check_combi(pepseq[p], mode) == target:
            result[allele].append(binder[allele][num][:, p])
            #new_array[:, w] = binder[allele][num][:, w]

        #result[allele].append(new_array)

    return result


def check_combi(pep, mode):
    target1 = aa_property[aa_property[mode] == 1]['aa'].tolist()
    target2 = aa_property[aa_property[mode] == 0]['aa'].tolist()

    if pep in target1:
        return 0
    elif pep in target2:
        return 1


def find_property_value(aa):
    try:
        return aa_property.loc[aa_property['aa'] == aa.upper()].values[0][2:5]
    except:
        return []

@ray.remote
def find_property2(df, target_group, binder, allele, target, mode, p):
    cp_value = {}
    cor_result = {}
    cor_result[allele] = []
    cp_value[allele] = []
    df = df[df['allele'].isin(target_group)]  # HLA-A,B,C 각각 가져오는 부분
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        #if check_combi(pepseq[p], mode) == target:
        cor_result[allele].append(binder[allele][num][:, p])
        cp_value[allele].append(find_property_value(pepseq[p]))

    return cp_value, cor_result


# if __name__ == "main":
print(sys.argv[1])
print(sys.argv[2])

object_store_memory = int(0.6 * ray.utils.get_system_memory() // 10 ** 9 * 10 ** 9)
try:
    ray.init(address='auto')
except:
    ray.init(dashboard_host='0.0.0.0',
             log_to_driver=False,
             _plasma_directory="/tmp",
             object_store_memory=object_store_memory
             )


p9_binder = load_gradcam_result()
p9_binder_id = ray.put(p9_binder)
df = load_pep_seq()
hla = load_short_hla()
hla_id = ray.put(hla)
df_id = ray.put(df)
del p9_binder, df, hla


aa_property = pd.read_excel('Amino_acid_property.xlsx')
aa_property['hydro'] = aa_property['Hydrophobicity'].map(lambda x: 1 if x >= 0 else 0)
aa_property['bulky'] = aa_property['Bulkiness'].map(lambda x: 1 if x >= np.median(aa_property['Bulkiness']) else 0)
aa_property['polar'] = aa_property['Polarity'].map(lambda x: 1 if x >= np.median(aa_property['Polarity']) else 0)

item = [[sys.argv[1]], [sys.argv[2]], [0, 1]]

for allele, mode, target in list(product(*item)):
    print(allele, mode, target)
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
        if sys.argv[3] == '0':
            result = ray.get([find_property.remote(df_id, total_g, p9_binder_id, allele, target, mode, p)
                              for allele in total_g])

            print('Saving Result')
            with open(f'/home/jaeung/Research/MHC/{allele}_{mode}_{target}_position_{p+1}_gradcam_result.pkl',
                      'wb') as f:
                pickle.dump(result, f)
        else:
            result = ray.get([find_property2.remote(df_id, total_g, p9_binder_id, allele, target, mode, p)
                              for allele in total_g])

            print('Saving Result')
            with open(f'/home/jaeung/Research/MHC/{allele}_{mode}_{target}_position_{p+1}_gradcam_result_with_cp_value.pkl',
                      'wb') as f:
                pickle.dump(result, f)
