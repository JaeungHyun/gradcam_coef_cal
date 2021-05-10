import pandas as pd
import numpy as np
import pickle
from itertools import combinations
from itertools import product
from tqdm import tqdm
import ray
import gc


@ray.remote
def cal_coef_by_p(binder1, binder2, p):
    rvalue = [np.corrcoef(binder1[num_i][:,p].reshape(-1),
                    binder2[num_j][:,p].reshape(-1))[0,1] \
              for num_i in range(len(binder1)) \
              for num_j in range(len(binder2))  ]
    
    rvalue = np.array(rvalue)
    rvalue = rvalue[~np.isnan(rvalue)]

    return rvalue


ray.init(dashboard_host='0.0.0.0')


with open('/home/jaeung/Research/MHC/ms+ba_short_hla_gradcam_result_remove_nan.pkl', 'rb') as f:
    p9_binder,p9_nonbinder, p10_binder, p10_nonbinder = pickle.load(f)


allele_list = list(p9_binder.keys())


group1 = ['HLA-A-2403', 'HLA-A-2402','HLA-A-2413','HLA-A-2301','HLA-A-2406','HLA-A-2407',]
group2 = ['HLA-A-3303', 'HLA-A-3301','HLA-A-6801','HLA-A-6601','HLA-A-3401','HLA-A-6602',
          'HLA-A-3101','HLA-A-7401',]
group3 = ['HLA-A-3001', 'HLA-A-0301', 'HLA-A-1101', 'HLA-A-1102','HLA-A-6812']
group4 = ['HLA-A-6802','HLA-A-6901']
group5 = ['HLA-A-0205','HLA-A-0206','HLA-A-0217','HLA-A-0216','HLA-A-0212','HLA-A-0219',
          'HLA-A-0207','HLA-A-0203','HLA-A-0201','HLA-A-0211','HLA-A-0204','HLA-A-0202',]
group6 = ['HLA-A-2601','HLA-A-2501','HLA-A-2608','HLA-A-2603','HLA-A-2602',]
group7 = ['HLA-A-0103','HLA-A-0101','HLA-A-2902','HLA-A-3002','HLA-A-3601','HLA-A-8001',]


for g in tqdm(['group1', 'group2', 'group3', 'group4', 'group5', 'group6', 'group7']):
    group1_ungroup = list(set( pd.Series(allele_list)[pd.Series(allele_list).str.contains('HLA-A')])\
                      -set(globals()[f'{g}']))
                     
    df_list = list(product(globals()[f'{g}'], group1_ungroup))
    for p in range(9):
        results = ray.get([cal_coef_by_p.remote(p9_binder[set1], p9_binder[set2], p) for set1, set2 in df_list])
    
        with open(f'/home/jaeung/Research/MHC/clustermap_correlation/short_HLA_A_{g}_P{p+1}_outgroup.pkl', 'wb') as f:
            pickle.dump(results, f)
        
        del results 
        gc.collect()