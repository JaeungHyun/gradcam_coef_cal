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
def find_property(df, target_group, binder, allele, target, mode):
    result = {}
    result[allele] = []
    df = df[df['allele'].isin(target_group)]
    for num, pepseq in enumerate(df.loc[df['allele'] == allele]['Peptide seq']):
        new_array = np.zeros((hla_len, 9))
        #for h, i in enumerate(hla[allele]):
        for w, j in enumerate(pepseq):
            if check_combi(j, mode) == target:
                new_array[:, w] = binder[allele][num][:, w]

        result[allele].append(new_array)

    return result


# if __name__ == "main":
print(sys.argv[1])
print(sys.argv[2])


try:
    ray.init(address='auto', log_to_driver=False)
except:
    ray.init(dashboard_host='0.0.0.0', log_to_driver=False)

if sys.argv[3]:
    p9_binder, _, _, _ = load_gradcam_result()
    p9_binder_id = ray.put(p9_binder)
    df = load_pep_seq()
    hla = load_short_hla()
    hla_id = ray.put(hla)
    df_id = ray.put(df)
    del p9_binder, df, hla

    with open('current_binder_id.pkl', 'wb') as f:
        pickle.dump((p9_binder_id, hla_id, df_id), f)
else:
    with open('current_binder_id.pkl', 'rb') as f:
        p9_binder_id, hla_id, df_id = pickle.load(f)


aa_property = pd.read_excel('Amino_acid_property.xlsx')
aa_property['hydro'] = aa_property['Hydrophobicity'].map(lambda x: 1 if x >= 0 else 0)
aa_property['bulky'] = aa_property['Bulkiness'].map(lambda x: 1 if x >= 15.4 else 0)  # 평균값이 15.367500000000001

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

    result = ray.get([find_property.remote(df_id, total_g, p9_binder_id, allele, target, mode)
                      for allele in total_g])

    print('Saving Result')
    with open(f'/home/jaeung/Research/MHC/{allele}_{mode}_{target}_gradcam_result.pkl', 'wb') as f:
        pickle.dump(result, f)
