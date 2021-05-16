import pandas as pd
import numpy as np
import pickle
from itertools import combinations
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

try:
    ray.init(dashboard_host='0.0.0.0',
             address='auto'
             )
except:
    ray.init(dashboard_host='0.0.0.0')

p9_binder, _, _, _ = load_gradcam_result()
allele_list = list(p9_binder.keys())
target_list, group_list = call_group_list(allele)
p9_binder_id = ray.put(p9_binder)


for g in tqdm(target_list):
    group_list = list(combinations(globals()[f'{g}'], 2))
    for p in range(9):
        results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])

        with open(f'/home/jaeung/Research/MHC/clustermap_correlation/short_HLA_A_{g}_P{p+1}_ingroup.pkl', 'wb') as f:
            pickle.dump(results, f)

        del results 
        gc.collect()
