import pandas as pd
import numpy as np
import re
import math
# import pysptools.distance.SAM
from pysptools.distance import SAM
import os
import warnings
warnings.filterwarnings('ignore')
import pickle

def get_SA(spec1 : list, spec2 : list)->float:
    mass_tolerance = 0.05
    spec1_length = len(spec1)
    spec2_length = len(spec2)
    
    #spectrum index
    idx1 = 0
    idx2 = 0
    
    #vector for mactched peak
    matched_vector1 = []
    matched_vector2 = []
    
    #spectrum size for normalization
    size1 = np.sqrt(sum([intensity[1]**2 for intensity in spec1]))
    size2 = np.sqrt(sum([intensity[1]**2 for intensity in spec2]))
    
    while(idx1 < spec1_length and idx2 < spec2_length):
        mz1, int1 = spec1[idx1]
        mz2, int2 = spec2[idx2]
        if abs(mz1 - mz2) <= mass_tolerance:
            matched_vector1.append(int1)
            matched_vector2.append(int2)
            idx1 += 1
            idx2 += 1
        elif mz1 > mz2:
            idx2 += 1
        else:
            idx1 += 1
    # size1 = np.sqrt(sum([intensity**2 for intensity in matched_vector1]))
    # size2 = np.sqrt(sum([intensity**2 for intensity in matched_vector2]))
    SA_val = 0
    
    try:
        SA_val = math.acos(np.dot(matched_vector1, matched_vector2)/(size1*size2))
    except:
        SA_val =  0
        
    SA_score = 1 - 2*(SA_val/np.pi)
    # print(f"SA matched pick length = {len(matched_vector1)}")
    return SA_score

def get_modified_SA(spec1 : list, spec2 : list, modified_mass, charge=1)->float:
    mass_tolerance = 0.05
    spec1_length = len(spec1)
    spec2_length = len(spec2)
    
    #spectrum index
    idx1 = 0
    idx2 = 0
    
    #vector for mactched peak
    matched_vector1 = []
    matched_vector2 = []
    
    #spectrum size for normalization
    size1 = np.sqrt(sum([intensity[1]**2 for intensity in spec1]))
    size2 = np.sqrt(sum([intensity[1]**2 for intensity in spec2]))
    # used_mz2_list = []
    while(idx1 < spec1_length and idx2 < spec2_length):
        
        mz1, int1 = spec1[idx1]
        mz2, int2 = spec2[idx2]
        if abs(mz1 - mz2) <= mass_tolerance:
            matched_vector1.append(int1)
            matched_vector2.append(int2)
            # used_mz2_list.append(mz2)
            idx1 += 1
            idx2 += 1
        elif mz1 > mz2:
            idx2 += 1
        else:
            idx1 += 1
    prev_vec1 = matched_vector1[:]
    prev_vec2 = matched_vector2[:]
    
    #find modified matching peak by mass shift
    idx1 = 0
    idx2 = 0
    while(idx1 < spec1_length and idx2 < spec2_length):
        mz1, int1 = spec1[idx1]
        if mz1 - modified_mass < 0:
            mz1 = 0
            # print("modifiedamss + mz1 i less than 0")
        else:
            mz1 -= modified_mass/charge
        mz2, int2 = spec2[idx2]
        if abs(mz1 - mz2) <= mass_tolerance and int1 not in matched_vector1:
            matched_vector1.append(int1)
            matched_vector2.append(int2)
            # used_mz2_list.append(mz2)
            idx1 += 1
            idx2 += 1
        elif abs(mz1-mz2) <= mass_tolerance and int1 in matched_vector1:
            #greedy하게 더 높은 픽을 가져감.
            if matched_vector2[matched_vector1.index(int1)] * int1 < int2*int1:
                matched_vector2[matched_vector1.index(int1)] = int2
                # print("greedy")
                
            idx1 += 1
            idx2 += 1
        elif mz1 > mz2:
            idx2 += 1
        else:
            idx1 += 1
            
    #find modified matching peak by mass shift charge = 2
    idx1 = 0
    idx2 = 0
    while(idx1 < spec1_length and idx2 < spec2_length):
        mz1, int1 = spec1[idx1]
        if mz1 - modified_mass < 0:
            mz1 = 0
            # print("modifiedamss + mz1 i less than 0")
        else:
            mz1 -= modified_mass/2
        mz2, int2 = spec2[idx2]
        if abs(mz1 - mz2) <= mass_tolerance and int1 not in matched_vector1:
            matched_vector1.append(int1)
            matched_vector2.append(int2)
            idx1 += 1
            idx2 += 1
        elif abs(mz1-mz2) <= mass_tolerance and int1 in matched_vector1:
            #greedy하게 더 높은 픽을 가져감.
            if matched_vector2[matched_vector1.index(int1)] * int1 < int2*int1:
                matched_vector2[matched_vector1.index(int1)] = int2
                # print("greedy")
                
            idx1 += 1
            idx2 += 1
        elif mz1 > mz2:
            idx2 += 1
        else:
            idx1 += 1
    SA_val = 0
    # size1 = np.sqrt(sum([intensity**2 for intensity in matched_vector1]))
    # size2 = np.sqrt(sum([intensity**2 for intensity in matched_vector2]))
    try:
        original_cosine = np.dot(prev_vec1, prev_vec2)/(size1*size2)
        cosine = np.dot(matched_vector1, matched_vector2)/(size1*size2)
        SA_val = math.acos(cosine)
        # print("done")
    except:
        SA_val  = 0
    SA_score = 1 - 2*(SA_val/np.pi)
    # print(f"modified SA matched pick length = {len(matched_vector1)}")
    return SA_score

# exit()
    
#SA를 같은 representative seqeunce의 모임에서 SA를 pair wise하게 구함
# grpup 1 unmodified seqeunce 끼리, group2 unmodified sequence끼리
# 그 이후에 SA들의 평균값을 취한뒤에 SA_mean_group1/2 컬럼 두개를 생성


# all_mgf = []
# with open("./MS2_embedding/data/SA/mgf.pkl", "rb") as f:
#     all_mgf = pickle.load(f)
# print(len(all_mgf[0][0]))
# print(len(all_mgf))


# exit()
# path = "./MS2_embedding/hek293_mgf/"
# file_list = os.listdir(path)
# print(file_list)
# all_mgf = []
# for mgf in file_list:
#     mgf_file = open(f"./MS2_embedding/hek293_mgf/{mgf}", "r")
#     contents = mgf_file.readlines()
#     mgf_file.close()
#     #indexing
#     index_begin = []
#     index_end = []
    
#     for i,line in enumerate(contents):
        
#         if line == "BEGIN IONS\n":
#             index_begin.append(i)
            
#         elif line == "END IONS\n":
#             index_end.append(i)     
#     print(mgf)       
#     all_mgf.append((contents, index_begin, index_end))

# print(len(all_mgf))
# with open("mgf.pkl", "wb") as f:
#     pickle.dump(all_mgf, f)
# exit()
# exit()
def get_spectrum(spec1_path, spec1_idx, mgf):
    spec1_path = spec1_path.split("\/")[-1]
    spec1_path = spec1_path.split("_")[4]
    path_int , path_letter = spec1_path[:2], spec1_path[2]
    
    if path_letter == 'B':
        path_int = int(path_int) + 11
    else:
        path_int = int(path_int) - 1
    
    spec1_mgf, index_begin, index_end = mgf[path_int]
    spec1_mgf = spec1_mgf[index_begin[spec1_idx-1]:index_end[spec1_idx-1]]
        
    spec1_list = []
    
    for a in spec1_mgf[5:]:
        spec1_list.append(list(map(float,a.split(" "))))
        
    return spec1_list
  
def get_pair_wise_SA():
    df = pd.read_csv("annotation_processed_input.csv")
    # print(df.shape)

    # df = df[:20000]
    # df = df.loc[df["Probability"] > 0.9, :]
    # print(df.shape)
    # df.to_csv("start_with_prob_filtered.csv", index=False)
    # exit()
    start_idx = 0
    prev_seq = df["seq_stripped"][0]
    result = []
    modification_list = []
    temp_list_for_unmod = []
    temp_list_for_mod = {}
    
    for i,(cur_seq, spec_path, file_idx, idx) in enumerate(zip(df["seq_stripped"],df["SpectrumFile"], df["Index"], df["idx"])):
        end_idx = i
        # print(prev_seq, cur_seq)
        if prev_seq != cur_seq or i == len(df):
            end_idx = i
            SA_list = []
            SA_mod_list = []
            SA_same_mod_list = []
            SA_same_mod_list_mod = []
            SA_rand_list = []
            dic = {}
            mod = 0
            
            for query_idx in range(start_idx, end_idx):
                #쿼리가 unmodified라면
                SA_mean = -1
                SA_mean2 = -1
                SA_mean3 = -1
                SA_mean4 = -1            
                if df.loc[query_idx]["ModificationAnnotation"] == '-1':
                    for cal_idx in range(start_idx, end_idx):
                        if df.loc[cal_idx]["ModificationAnnotation"] == '-1':
                            
                            if cal_idx == query_idx:
                                continue
                            
                            spec1 = get_spectrum(df.loc[query_idx]["SpectrumFile"], df.loc[query_idx]["Index"],all_mgf)
                            spec2 = get_spectrum(df.loc[cal_idx]["SpectrumFile"], df.loc[cal_idx]["Index"],all_mgf)
                            SA_list.append(get_SA(spec1, spec2))
                            temp_list_for_unmod.append(get_SA(spec1, spec2))
                    if SA_list:
                        SA_mean = max(SA_list)
                    else:
                        SA_mean = None
                    mod = 0
                          
                #query가 modified라면
                else:
                    for cal_idx in range(start_idx, end_idx):
                        #cal peptide 가 modified가 아닐때 : mod - umod 비교
                        if df.loc[cal_idx]["ModificationAnnotation"] == '-1':
                            if cal_idx == query_idx:
                                continue
                            
                            # if i > 600:
                            #     print(df.loc[query_idx]["Peptide"])
                            #     print(df.loc[cal_idx]["Peptide"])
                            spec1 = get_spectrum(df.loc[query_idx]["SpectrumFile"], df.loc[query_idx]["Index"],all_mgf)
                            spec2 = get_spectrum(df.loc[cal_idx]["SpectrumFile"], df.loc[cal_idx]["Index"],all_mgf)
                            # print("---------------------start")
                            delta_mass = df.loc[query_idx]["CalculatedMW"] - df.loc[cal_idx]["CalculatedMW"]        
                            # if "C+" in df.loc[query_idx]["Peptide"]:
                            #     if get_modified_SA(spec1, spec2, df.loc[query_idx]["modified_mass"]-57.021) >get_modified_SA(spec1, spec2, df.loc[query_idx]["modified_mass"]):
                            #         SA_mod_list.append(get_modified_SA(spec1, spec2, df.loc[query_idx]["modified_mass"]-57.021))
                            #     else:
                            #         SA_mod_list.append(get_modified_SA(spec1, spec2, df.loc[query_idx]["modified_mass"]))
                            # else:
                            if df.loc[query_idx]["mod_processed"] not in temp_list_for_mod:
                                temp_list_for_mod[df.loc[query_idx]["mod_processed"]] = [get_modified_SA(spec1, spec2, delta_mass)]
                            else:
                                temp_list_for_mod[df.loc[query_idx]["mod_processed"]].append(get_modified_SA(spec1, spec2, delta_mass))
                            # temp_list_for_mod.append((get_modified_SA(spec1, spec2, delta_mass), df.loc[query_idx]["ModificationAnnotation"]))
                            SA_mod_list.append(get_modified_SA(spec1, spec2, delta_mass)) #
                            # SA_list.append(df.loc[query_idx]["Peptide"])
                            SA_list.append(get_SA(spec1, spec2))
                            # print("---------------------end")
                        # mod-mod 비교
                        else:
                            if cal_idx == query_idx:
                                continue
                            if df.loc[query_idx]["Peptide"] == df.loc[cal_idx]["Peptide"]:
                                delta_mass = df.loc[query_idx]["CalculatedMW"] - df.loc[cal_idx]["CalculatedMW"]
                                spec1 = get_spectrum(df.loc[query_idx]["SpectrumFile"], df.loc[query_idx]["Index"],all_mgf)
                                spec2 = get_spectrum(df.loc[cal_idx]["SpectrumFile"], df.loc[cal_idx]["Index"],all_mgf)
                                # if df.loc[query_idx]["ModificationAnnotation"] is not np.nan:
                                    
                                # print("---------------------start with mod")
                                SA_same_mod_list_mod.append(get_modified_SA(spec1, spec2, delta_mass)) #
                                # SA_list.append(df.loc[query_idx]["Peptide"])
                                SA_same_mod_list.append(get_SA(spec1, spec2)) 
                                # print("---------------------end with mod")
                            
                    if SA_list:
                        SA_mean = max(SA_list)
                    else:
                        SA_mean = None
                    if SA_mod_list:
                        SA_mean2 = max(SA_mod_list)
                    else:
                        SA_mean2 = None
                    if SA_same_mod_list:
                        SA_mean3 = max(SA_same_mod_list)
                    else:
                        SA_mean3 = None
                    if SA_same_mod_list_mod:
                        SA_mean4 = max(SA_same_mod_list_mod)
                    else:
                        SA_mean4 = None
                    mod = 1
                SA_list = []
                SA_mod_list = []
                SA_same_mod_list = []
                SA_same_mod_list_mod = []
                
                #query sequence에 대한 random match와의 SA
                query_mass = df.loc[query_idx]["ObservedMW"]
                
                df_new = df.loc[abs(df["ObservedMW"] - query_mass) < 350,["idx"]]
                rand_idx = df_new.sample(n=1)["idx"].item()
                
                spec1 = get_spectrum(df.loc[query_idx]["SpectrumFile"], df.loc[query_idx]["Index"],all_mgf)
                series= df.loc[df["idx"] == rand_idx,["SpectrumFile", "Index", "CalculatedMW"]]
                spec2_file = series["SpectrumFile"].item()
                spec2_idx = series["Index"].item()
                
                spec2 = get_spectrum(spec2_file, spec2_idx,all_mgf)
                # SA_val = get_SA(spec1, spec2)
                SA_val = get_modified_SA(spec1, spec2, df.loc[query_idx]["CalculatedMW"] - series["CalculatedMW"].item())
                result.append((mod,SA_mean, SA_mean2, SA_mean3, SA_mean4, df.loc[query_idx]["idx"], 2, SA_val, rand_idx,start_idx,df.loc[query_idx]["ModificationAnnotation"], df.loc[query_idx]["SpectrumFile"], df.loc[query_idx]["Index"], df.loc[query_idx]["seq_stripped"],df.loc[query_idx]["CalculatedMW"])) 
            
            for mod_annotation in temp_list_for_mod:
                if temp_list_for_unmod:
                    modification_list.append((df.loc[start_idx]["seq_stripped"], mod_annotation, max(temp_list_for_unmod), max(temp_list_for_mod[mod_annotation])))
            prev_seq = cur_seq
            start_idx = i
            temp_list_for_unmod = []
            temp_list_for_mod = {}
        # if len(result) > 5:
        #     print(i, result[-4:])
        if i % 1000 == 0:
            print(i)
    with open("modification_class","w") as f:
        f.write("seq_stripped,Annotation,unmod-unmod,mod-unmod\n")
        for a in modification_list:
            f.write(f"{a[0]},{a[1]},{a[2]},{a[3]}\n")
    with open("SA_mean.txt", "w") as f:
        f.write("mod,SA,SA_mod,SA_same_mod,SA_same_mod_with_mod,idx,rand,rand_SA,rand_idx,start_idx,ModificationAnnotation,SpectrumFile,Index,seq_stripped,CalculatedMW\n")
        for i,a in enumerate(result):
            f.write(f"{a[0]},{a[1]},{a[2]},{a[3]},{a[4]},{a[5]},{a[6]},{a[7]},{a[8]},{a[9]},{a[10]},{a[11]},{a[12]},{a[13]},{a[14]}\n")
                   
        
# get_pair_wise_SA()
