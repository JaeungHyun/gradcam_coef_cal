import numpy as np
import sys
from util import *
import pickle


allele = sys.argv[1]
mode = sys.argv[2]  #뭐였더라..
group_mode = sys.argv[3]  # ingroup or outgroup

for p in range(9):
    print('importing binder data')
    data_list = load_target_gradcam_result(allele, mode, 0, p, cp=None)  # 어짜피 polar안에 다 있음 다른 cp들 결과

    if mode == 'hydro':
        v = 0
    elif mode == 'bulky':
        v = 1
    elif mode == 'polar':
        v = 2
    elif mode == 'Charge':
        v = 3
    elif mode == 'MW':
        v = 4

    key_list = []
    for data in data_list:
        #key_list.extend(list(data.keys()))
        for list_ in data:
            key_list.extend(list(list_.keys()))

    cp_result = {}
    result = {}
    for key in key_list:
        cp_result[key] = []
        result[key] = []

    print('cp result')
    for data in data_list:
        #key_lists = []
        #for key in data[0].keys():
        for key in data.keys():
            #for i, value in enumerate(data[key]):
            for i, value in enumerate(data[0][key]):
                try:
                    cp_result[key].append(value)
                except:
                    print(value)
        for key in data[1].keys():
            for value in data[1][key]:
                result[key].append(value)
    with open(f'/home/jaeung/Research/MHC/short_{allele}_{mode}_{group_mode}_{p+1}_with_gradcam_result.pkl', 'wb') as f:
        pickle.dump((cp_result, result), f)
        #pickle.dump((cp_result), f)
