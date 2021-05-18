import pandas as pd
import numpy as np
import pickle
from itertools import combinations, product
from tqdm import tqdm
import ray
import sys
from util import *


@ray.remote
def cal_coef_by_p(binder_id, binder1, binder2, p):
    rvalue = [np.corrcoef(binder_id[binder1][num_i][:,p].reshape(-1),
                          binder_id[binder2][num_j][:,p].reshape(-1))[0,1] \
              for num_i in range(len(binder_id[binder1])) \
              for num_j in range(len(binder_id[binder2]))  ]
    
    rvalue = np.array(rvalue)
    rvalue = rvalue[~np.isnan(rvalue)]

    return rvalue

'''
    Ray를 터미널에서 실행하고 스크립트 실행할시 
    address='auto'로 하면 자동으로 이미 생성되어있는 곳에 들어감
    클러스터가 여러개면 ip를 지정해주어야함
'''
allele = sys.argv[1]
mode = sys.argv[2]  # total, [hydro, bulky, 0,1,2,3] 나누기, 2021.05.16 현재 hydro, bulky 0~3만 계산하면됨
group_mode = sys.argv[3]  # ingroup or outgroup


try:
    ray.init(dashboard_host='0.0.0.0',
             address='auto'
             )
except:
    ray.init(dashboard_host='0.0.0.0')

target_list, target_group_list = call_group_list(allele)

for target in range(4):
    p9_binder = load_target_gradcam_result(allele, mode, target)
    if mode != 'total':
        result = {}
        for dic in p9_binder:
            for key, value in dic.items():
                result[key] = value
    p9_binder = result

    p9_binder_id = ray.put(p9_binder)
    allele_list = list(p9_binder.keys())

    for i, g in tqdm(enumerate(target_list)):
        if group_mode == 'ingroup':
            group_list = list(combinations(target_group_list[i], 2))

        elif group_mode == 'outgroup':
            outgroup = pd.Series(allele_list)[pd.Series(allele_list).str.contains(f'{allele}')] - target_group_list[i]
            group_list = list(product(target_group_list[i], outgroup))

        for p in range(9):
            results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])
            with open(f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{target}_{g}_P{p+1}_{group_mode}.pkl', 'wb') as f:
                pickle.dump(results, f)

            del results
            gc.collect()
