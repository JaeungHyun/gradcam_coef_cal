import numpy as np
import sys
from util import *
import pickle


allele = sys.argv[1]
mode = sys.argv[2]  # total, [hydro, bulky, 0,1,2,3] 나누기, 2021.05.16 현재 hydro, bulky 0~3만 계산하면됨
group_mode = sys.argv[3]  # ingroup or outgroup



for p in range(9):
    print('importing binder data')
    data_list = load_target_gradcam_result(allele, 'polar', 0, p, cp='cp')  # 어짜피 polar안에 다 있음 다른 cp들 결과
    if mode == 'hydro':
        v = 0
    elif mode == 'bulky':
        v = 1
    elif mode == 'polar':
        v = 2

    key_list = []
    for data in data_list:
        for list_ in data:
            key_list.extend(list(list_.keys()))

    #print(key_list)

    cp_result = {}
    result = {}
    for key in key_list:
        cp_result[key] = []
        result[key] = []
    #print(cp_result)
    print('cp result')
    for data in data_list:
        key_lists = []
        for key in data[0].keys():
            for i, value in enumerate(data[0][key]):
                try:
                    # print(type(value[v]))
                #print(key)

                    cp_result[key].append(value[v])
                ##key_lists.append(key)
                except:
                    print(value)
                #    print(key, i)

        for key in data[1].keys():
            #print(data[1].keys())
            for value in data[1][key]:
                result[key].append(value.astype('float16'))
    with open(f'/home/jaeung/Research/MHC/short_{allele}_{mode}_{group_mode}_{p+1}_with_gradcam_result.pkl',
                'wb') as f:
            pickle.dump((cp_result, result), f)