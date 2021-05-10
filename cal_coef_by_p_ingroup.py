import pandas as pd
import numpy as np
import pickle
from itertools import combinations
from itertools import product
from tqdm import tqdm
import ray
import gc


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
ray.init(dashboard_host='0.0.0.0',
         #address='163.180.14.45:6379'
         )


with open('/home/jaeung/Research/MHC/ms+ba_short_hla_gradcam_result_remove_nan.pkl', 'rb') as f:
    p9_binder, _, _, _ = pickle.load(f)


allele_list = list(p9_binder.keys())


p9_binder_id = ray.put(p9_binder)

allele = 'HLA-C'

if allele == 'HLA-A':
    group1 = ['HLA-A-2403', 'HLA-A-2402', 'HLA-A-2413', 'HLA-A-2301', 'HLA-A-2406', 'HLA-A-2407']
    group2 = ['HLA-A-3303', 'HLA-A-3301', 'HLA-A-6801', 'HLA-A-6601', 'HLA-A-3401', 'HLA-A-6602',
              'HLA-A-3101', 'HLA-A-7401']
    group3 = ['HLA-A-3001', 'HLA-A-0301', 'HLA-A-1101', 'HLA-A-1102', 'HLA-A-6812']
    group4 = ['HLA-A-6802', 'HLA-A-6901']
    group5 = ['HLA-A-0205', 'HLA-A-0206', 'HLA-A-0217', 'HLA-A-0216', 'HLA-A-0212', 'HLA-A-0219',
              'HLA-A-0207', 'HLA-A-0203', 'HLA-A-0201', 'HLA-A-0211', 'HLA-A-0204', 'HLA-A-0202']
    group6 = ['HLA-A-2601', 'HLA-A-2501', 'HLA-A-2608', 'HLA-A-2603', 'HLA-A-2602']
    group7 = ['HLA-A-0103', 'HLA-A-0101', 'HLA-A-2902', 'HLA-A-3002', 'HLA-A-3601', 'HLA-A-8001']

    target_list = ['group1', 'group2', 'group3', 'group4', 'group5', 'group6', 'group7']

elif allele == 'HLA-B':
    group1 = ['HLA-B-5301', 'HLA-B-3501', 'HLA-B-3507', 'HLA-B-3508', 'HLA-B-1511']
    group2 = ['HLA-B-0704', 'HLA-B-0702', 'HLA-B-4201', 'HLA-B-3502', 'HLA-B-3503','HLA-B-3504','HLA-B-3506',]
    group3 = ['HLA-B-8101', 'HLA-B-4202',]
    group4 = ['HLA-B-5401','HLA-B-5501', 'HLA-B-5502','HLA-B-5601']
    group5 = ['HLA-B-5101', 'HLA-B-5108', 'HLA-B-7301', 'HLA-B-3906',]
    group6 = ['HLA-B-2710', 'HLA-B-2702', 'HLA-B-2701', 'HLA-B-2704', 'HLA-B-2703', 'HLA-B-2705', 'HLA-B-2708',
              'HLA-B-2707', 'HLA-B-2706',]
    group7 = ['HLA-B-3905', 'HLA-B-3901', 'HLA-B-3801', 'HLA-B-3802', 'HLA-B-1509', 'HLA-B-1510',]
    group8 = ['HLA-B-3924', 'HLA-B-1402', 'HLA-B-1403',]
    group9 = ['HLA-B-2709','HLA-B-3909',]
    group10 = ['HLA-B-4901', 'HLA-B-5001', 'HLA-B-4006', 'HLA-B-4101', 'HLA-B-4501',]
    group11 = ['HLA-B-1803', 'HLA-B-1801', 'HLA-B-4402', 'HLA-B-4403', 'HLA-B-4427', 'HLA-B-4428',]
    group12 = ['HLA-B-4102', 'HLA-B-4104', 'HLA-B-4103', 'HLA-B-4409', 'HLA-B-4002', 'HLA-B-4001',]
    group13 = ['HLA-B-1508','HLA-B-1501','HLA-B-1503','HLA-B-1502','HLA-B-4601',]
    group14 = ['HLA-B-5703','HLA-B-5701','HLA-B-5801', 'HLA-B-5802','HLA-B-1517',]
    group15 = ['HLA-B-5201', 'HLA-B-1302', ]
    group16 = ['HLA-B-0803', 'HLA-B-0802',]

    target_list = ['group1', 'group2', 'group3', 'group4', 'group5',
                   'group6', 'group7', 'group8', 'group9', 'group10',
                   'group11', 'group12', 'group13', 'group14', 'group15',
                   'group16']

else:
    group1 = ['HLA-C-0401', 'HLA-C-0501', 'HLA-C-0403', 'HLA-C-0802', ]
    group2 = ['HLA-C-1402', 'HLA-C-1402', ]
    group3 = ['HLA-C-0704', 'HLA-C-0702', 'HLA-C-0602', 'HLA-C-0701', ]
    group4 = ['HLA-C-1502', 'HLA-C-1505', ]
    group5 = ['HLA-C-1701', 'HLA-C-0801', 'HLA-C-0304', 'HLA-C-0303', 'HLA-C-1202',
              'HLA-C-0202', 'HLA-C-1203', 'HLA-C-1601', 'HLA-C-0302', ]
    target_list = ['group1', 'group2', 'group3', 'group4', 'group5']


for g in tqdm(target_list):
    group_list = list(combinations(globals()[f'{g}'], 2))
    for p in range(9):
        results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])

        with open(f'/home/jaeung/Research/MHC/clustermap_correlation/short_HLA_A_{g}_P{p+1}_ingroup.pkl', 'wb') as f:
            pickle.dump(results, f)

        del results 
        gc.collect()
