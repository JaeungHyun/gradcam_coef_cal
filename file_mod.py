import numpy as np
import sys
from util import *
import pickle
from itertools import product


item = [[sys.argv[1]], ['polar','hydro','bulky', 'MW', 'Charge'], ['random', 'natural']]

for allele, mode, false_Kinds in list(product(*item)):
    for p in range(9):
        print('importing binder data')
        data_list = load_target_gradcam_result(allele, mode, p, false_Kinds)

        key_list = []
        if mode == 'polar':
            for data in data_list:
                for list_ in data:
                    key_list.extend(list(list_.keys()))

            cp_result = {}
            result = {}
            for key in key_list:
                cp_result[key] = []
                result[key] = []
            print('cp result')
            for data in data_list:
                for key in data[0].keys():
                    for i, value in enumerate(data[0][key]):
                        try:
                            cp_result[key].append(value)
                        except:
                            print(value)
                for key in data[1].keys():
                    for value in data[1][key]:
                        result[key].append(value)

            with open(f'/home/jaeung/960evo/result/by_position/short_{allele}_{p+1}_with_gradcam_by_position.pkl', 'wb') as f:
                pickle.dump(cp_result, f)
            with open(f'/home/jaeung/960evo/result/by_position/short_{allele}_{mode}_{p+1}_with_cp_value.pkl', 'wb') as f:
                pickle.dump(result, f)
        else:
            cp_result = {}
            for data in data_list:
                key_list.extend(list(data.keys()))
            for key in key_list:
                cp_result[key] = []

            print('cp result')
            for data in data_list:

                for key in data.keys():
                    for i, value in enumerate(data[key]):
                        try:
                            cp_result[key].append(value)
                        except:
                            print(value)
            with open(f'/home/jaeung/960evo/result/short_{allele}_{mode}_{p + 1}_with_cp_value.pkl', 'wb') as f:
                pickle.dump(cp_result, f)


