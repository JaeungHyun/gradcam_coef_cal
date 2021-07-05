import numpy as np
from tqdm import tqdm
import ray
import sys
from util import *


@ray.remote
def cal_coef_by_p(binder_id, binder1, binder2):
    rvalue = [np.corrcoef(binder_id[binder1][num_i].reshape(-1),
                          binder_id[binder2][num_j].reshape(-1))[0, 1]\
              for num_i in range(len(binder_id[binder1]))\
              for num_j in range(len(binder_id[binder2]))]
    
    rvalue = np.array(rvalue)
    rvalue = rvalue[~np.isnan(rvalue)]
    return rvalue


@ray.remote
def cal_coef_by_matrix(binder_id, binder1, binder2):
    rvalue = [np.corrcoef(binder_id[binder1][num_i].reshape(-1),
                          binder_id[binder2][num_j].reshape(-1))[0, 1]\
              for num_i in range(len(binder_id[binder1]))\
              for num_j in range(len(binder_id[binder2]))]

    rvalue = np.array(rvalue)
    rvalue = rvalue[~np.isnan(rvalue)]
    return rvalue


@ray.remote
def cal_coef_by_p_with_cp_value(cp_id, binder_id, allele1, allele2):
    sum_cp = [cp_id[allele1][num_i] + cp_id[allele2][num_j]
              for num_i in range(len(cp_id[allele1])) \
              for num_j in range(len(cp_id[allele2]))]
    mul_cp = [cp_id[allele1][num_i] * cp_id[allele2][num_j]
              for num_i in range(len(cp_id[allele1])) \
              for num_j in range(len(cp_id[allele2]))]
    rvalue = [np.corrcoef(binder_id[allele1][num_i],
                          binder_id[allele2][num_j])[0, 1] \
              for num_i in range(len(binder_id[allele1])) \
              for num_j in range(len(binder_id[allele2]))]
    sum_cp = np.array(sum_cp)
    mul_cp = np.array(mul_cp)
    rvalue = np.array(rvalue)

    return sum_cp, mul_cp, rvalue

'''
    Ray를 터미널에서 실행하고 스크립트 실행할시 
    address='auto'로 하면 자동으로 이미 생성되어있는 곳에 들어감
    클러스터가 여러개면 ip를 지정해주어야함
'''
allele = sys.argv[1]
mode = sys.argv[2]  # total, [hydro, bulky, 0,1,2,3] 나누기, 2021.05.16 현재 hydro, bulky 0~3만 계산하면됨
group_mode = sys.argv[3]  # ingroup or outgroup
initial = sys.argv[4]

object_store_memory = int(0.6 * ray.utils.get_system_memory() // 10 ** 9 * 10 ** 9)


try:
    ray.init(address='auto', log_to_driver=False)
except:
    ray.init(dashboard_host='0.0.0.0',
             log_to_driver=False,
             _plasma_directory="/tmp",
             object_store_memory=object_store_memory
             )
target_list, target_group_list = call_group_list(allele)

if sys.argv[4] == "cp":
    for p in range(9):
        print('importing binder data')
        data_list = load_target_gradcam_result(allele, 'polar', 0, p, cp='cp')  # 어짜피 polar안에 다 있음 다른 cp들 결과
        # p9_binder = load_target_gradcam_result(allele, mode, target, p)
        cp_result = {}
        result = {}
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

        for key in key_list:
            cp_result[key] = []
            result[key] = []

        for data in data_list:
            key_list = []
            for key in data[0].keys():
                for i, value in enumerate(data[0][key]):
                    try:
                        cp_result[key].append(value[v])
                        key_list.append(key)
                    except:
                        print(key, i)
            for key in key_list:
                for value in data[1][key]:
                    result[key].append(value)

        cp_value_id = ray.put(cp_result)
        p9_binder_id = ray.put(result)
        allele_list = list(result.keys())
        del result, cp_result

        for i, g in tqdm(enumerate(target_list)):
            group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
            print(allele, mode, g, f'P{p + 1}\n')
            results = ray.get(
                [cal_coef_by_p_with_cp_value.remote(cp_value_id, p9_binder_id, set1, set2) for set1, set2 in
                 group_list])
            with open(
                    f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{g}_{group_mode}_{p+1}_with_cp_value.pkl',
                    'wb') as f:
                pickle.dump(results, f)

        del results

elif sys.argv[4] != "cp" and (mode == 'hydro' or mode == "bulky" or mode == 'polar'):
    for p in range(9):
        for target in range(2):
            p9_binder = load_target_gradcam_result(allele, mode, target, p)
            result = {}
            for dic in p9_binder:
                for key, value in dic.items():
                    result[key] = value

            p9_binder_id = ray.put(result)
            allele_list = list(result.keys())
            del p9_binder

            for i, g in tqdm(enumerate(target_list)):
                group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
                print(allele, mode, target, g, f'P{p + 1}\n')
                results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2) for set1, set2 in group_list])
                with open(
                        f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{target}_{g}_P{p + 1}_{group_mode}.pkl',
                        'wb') as f:
                    pickle.dump(results, f)

                del results
                gc.collect()
        del p9_binder_id

elif mode == "pattern":
    p9_binder = load_target_gradcam_result(allele, mode)
    print('importing binder data')
    p9_binder_id = ray.put(p9_binder)
    allele_list = list(p9_binder.keys())
    del p9_binder
    with open('total_binder_id.pkl', 'wb') as f:
        pickle.dump((p9_binder_id, allele_list), f)

    for i, g in tqdm(enumerate(target_list[10:13])):
        group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
        print(allele, mode, g)
        results = ray.get([cal_coef_by_matrix.remote(p9_binder_id, set1, set2) for set1, set2 in group_list])
        with open(f'short_{allele}_{mode}_{g}_{group_mode}.pkl', 'wb') as f:
            pickle.dump(results, f)

        del results
        gc.collect()


