import numpy as np
import pandas as pd
from tqdm import tqdm
import ray
import sys
from util import *
import pickle


@ray.remote
def cal_coef_by_matrix(binder_id, binder1, binder2):
    rvalue = [np.corrcoef(array1.reshape(-1),
                          array2.reshape(-1))[0, 1] \
              for array1 in binder_id[binder1] \
              for array2 in binder_id[binder2]]

    rvalue = np.array(rvalue)
    # rvalue = rvalue[~np.isnan(rvalue)]
    return rvalue


allele = sys.argv[1]
group_mode = sys.argv[2]  # ingroup or outgroup


try:
    ray.init(address='auto', log_to_driver=False)
except:
    ray.init(dashboard_host='0.0.0.0',
             log_to_driver=False,
             _plasma_directory="/tmp",
             )


target_list, target_group_list = call_group_list(allele)

item = [[sys.argv[1]], ['ingroup', 'outgroup'],
        ['random', 'natural']]

for allele, group_mode, false_kinds in list(product(*item)):
    p9_binder = load_gradcam_result(false_kinds)
    p9_binder_id = ray.put(p9_binder)
    allele_list = pd.Series(p9_binder.keys())
    allele_list = allele_list[allele_list.str.contains(allele)]

    del p9_binder

    target_list, group_list = call_group_list(allele)
    total_g = []
    for g in group_list:
        total_g.extend(g)
    if group_mode == 'ingroup':
        print('Calculating ingroup')
        for i, g in enumerate(tqdm(target_list)):
            group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
            results = ray.get(
                [cal_coef_by_matrix.remote(p9_binder_id, set1, set2) for set1, set2 in
                 group_list])

            np.save(f'/home/hdd/result/clustermap_correlation/{allele}_{g}_{group_mode}_{false_kinds}_gradcam_cor.pkl',
                    results, allow_pickle=True)
    else:
        print('Calculating outgroup')
        for i, g in enumerate(tqdm(target_list)):
            group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
            results = ray.get(
                [cal_coef_by_matrix.remote(p9_binder_id, set1, set2) for set1, set2 in
                 group_list])

            np.save(
                f'/home/hdd/result/clustermap_correlation/{allele}_{g}_{group_mode}_{false_kinds}_gradcam_cor.pkl',
                results, allow_pickle=True)
