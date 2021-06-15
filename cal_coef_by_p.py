import numpy as np
from tqdm import tqdm
import ray
import sys
from util import *


@ray.remote
def cal_coef_by_p(binder_id, binder1, binder2, p):
    rvalue = [np.corrcoef(binder_id[binder1][num_i][:,p].reshape(-1),
                          binder_id[binder2][num_j][:,p].reshape(-1))[0, 1]\
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


if mode == 'hydro' or mode == "bulky":
    for target in range(2):
        p9_binder = load_target_gradcam_result(allele, mode, target)

        result = {}
        for dic in p9_binder:
            for key, value in dic.items():
                result[key] = value

        p9_binder_id = ray.put(result)
        allele_list = list(result.keys())
        del p9_binder

        for i, g in tqdm(enumerate(target_list)):
            group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)

            for p in range(9):
                print(allele, mode, target, g, f'P{p+1}')
                results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])
                with open(f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{target}_{g}_P{p+1}_{group_mode}.pkl', 'wb') as f:
                    pickle.dump(results, f)

                del results
                gc.collect()
        del p9_binder_id

elif mode == "pattern":
    if initial == '1':
        p9_binder = load_target_gradcam_result(allele, mode)
        print('importing binder data')
        p9_binder_id = ray.put(p9_binder)
        allele_list = list(p9_binder.keys())
        del p9_binder
        with open('total_binder_id.pkl', 'wb') as f:
            pickle.dump((p9_binder_id, allele_list), f)
    else:
        print('Use exist data')
        with open('total_binder_id.pkl', 'rb') as f:
            p9_binder_id, allele_list = pickle.load(f)

    for i, g in tqdm(enumerate(target_list[10:13])):
        group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
        print(allele, mode, g)
        results = ray.get([cal_coef_by_matrix.remote(p9_binder_id, set1, set2) for set1, set2 in group_list])
        with open(f'short_{allele}_{mode}_{g}_{group_mode}.pkl', 'wb') as f:
            pickle.dump(results, f)

        del results
        gc.collect()

else:
    p9_binder = load_target_gradcam_result(allele, mode)

    if initial == '1':
        print('importing binder data')
        p9_binder_id = ray.put(p9_binder)
        allele_list = list(p9_binder.keys())
        del p9_binder
        with open('total_binder_id.pkl', 'wb') as f:
            pickle.dump((p9_binder_id, allele_list), f)
    else:
        print('Use exist data')
        with open('total_binder_id.pkl', 'rb') as f:
            p9_binder_id, allele_list = pickle.load(f)

    for i, g in tqdm(enumerate(target_list)):
        group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)

        for p in range(9):
            print(allele, mode, g, f'P{p + 1}')
            results = ray.get([cal_coef_by_p.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])
            with open(
                f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}__{g}_P{p+1}_{group_mode}.pkl',
                    'wb') as f:
                pickle.dump(results, f)

            del results
            gc.collect()
    del p9_binder_id

