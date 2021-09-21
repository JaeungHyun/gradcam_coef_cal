import numpy as np
from tqdm import tqdm
import ray
import sys
from util import *
import pickle


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
def cal_coef_by_p_with_cp_value(binder_id, allele1, allele2):
    rvalue = [np.corrcoef(num_i, num_j)[0, 1] \
              for num_i in binder_id[allele1] \
              for num_j in binder_id[allele2]]
    rvalue = np.array(rvalue)

    return rvalue


@ray.remote
def cal_coef_by_p_with_cp_sub_value(cp_id, allele1, allele2):
    sub_cp = [1/(np.abs(cp_id[allele1][num_i] - cp_id[allele2][num_j]) + 1)
              for num_i in range(len(cp_id[allele1])) \
              for num_j in range(len(cp_id[allele2]))]

    sub_cp = np.array(sub_cp)
    return sub_cp


'''
    Ray를 터미널에서 실행하고6 스크립트 실행할시 
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
    if mode == 'polar':
        v=2
    elif mode == 'bulky':
        v=1
    for p in range(9):
        print('importing binder data')
        data = load_target_gradcam_result(allele, mode, 0, p, cp='cp') # cp 데이터 불러옴
        p9_binder = load_gradcam_result()
        #p9_binder_id = ray.put(p9_binder)

        cp_result = {}
        for key, value in data.items():
            cp_result[key] = value

        cp_value_id = ray.put(cp_result)
        allele_list = list(p9_binder.keys())
        del p9_binder, cp_result, data
        #del p9_binder

        for i, g in enumerate(tqdm(target_list)):
            group_list = return_group_list(group_mode, target_group_list, allele_list, allele, i)
            print(allele, mode, g, f'P{p + 1}\n')
            # print('CP sum')
            # results = ray.get(
            #     [cal_coef_by_p_with_cp_sum_value.remote(cp_value_id, set1, set2) for set1, set2 in
            #      group_list])
            # with open(
            #         f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{g}_{group_mode}_{p + 1}_with_cp_sum_value.pkl',
            #         'wb') as f:
            #     pickle.dump(results, f)
            # del results
            # print('abs Cor coef')
            # cor_results = ray.get([cal_coef_by_p_abs.remote(p9_binder_id, set1, set2, p) for set1, set2 in group_list])
            # with open(
            #         f'/home/jaeung/Research/MHC/short_{allele}_{mode}_{g}_{group_mode}_{p + 1}_with_cor_abs.pkl',
            #         'wb') as f:
            #     pickle.dump(cor_results, f)
            print('GradCAM correlation')
            results = ray.get(
                [cal_coef_by_p_with_cp_value.remote(cp_value_id, set1, set2) for set1, set2 in
                 group_list])
            with open(
                    f'/data/result/clustermap_correlation/short_{allele}_{mode}_{g}_{group_mode}_{p+1}_gradcam_cor.pkl',
                    'wb') as f:
                pickle.dump(results, f)
            del results

            print('CP subs')
            results = ray.get(
                [cal_coef_by_p_with_cp_sub_value.remote(cp_value_id, set1, set2) for set1, set2 in
                 group_list])
            with open(
                    f'/home/jaeung/Research/MHC/clustermap_correlation/short_{allele}_{mode}_{g}_{group_mode}_{p+1}_with_cp_sub_value.pkl',
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
                        f'/home/hdd/Research/clustermap_correlation/short_{allele}_{mode}_{target}_{g}_P{p + 1}_{group_mode}.pkl',
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


