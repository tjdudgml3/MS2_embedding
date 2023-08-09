import pandas as pd
import numpy as np
import re
# import pysptools.distance.SAM
from pysptools.distance import SAM


# modified sequence와 unmodified sequence를 그루핑하는 함수
def group_seqeunce():
    df = pd.read_csv("chick_modplus_fdr1.tsv", delimiter="\t")
    print(df.shape)
    df = df.loc[:,["SpectrumFile", "Index", "ObservedMW", "Peptide"]]
    print(df.shape)
    
    
    sequence_dic = {}
    #이중포문으로 grouping할것 실행시간 한시간반예상
    for idx1, query_seq in enumerate(df["Peptide"]):
        # modifided_mass = re.findall("[0-9]+\.*[0-9]*", query_seq)
        if idx1 % 10000 == 0:
            print(idx1)
        
        stripped_seq = re.sub("\{","", query_seq)
        stripped_seq = re.sub("\}","", stripped_seq)
        stripped_seq = re.sub("\+[0-9]+\.*[0-9]*","", stripped_seq)
        # stripped_seq = query_seq.replace(f"{}.*", "")
        # if idx1 > 100:
        #     print(sequence_dic)
        #     break
        # modified_mass_sum =  sum(map(float,modifided_mass))
        # print(modifided_mass, modified_mass_sum)
        # print(stripped_seq, query_seq)
        if stripped_seq not in sequence_dic:
            sequence_dic[stripped_seq] = []#(query_seq, modified_mass_sum, idx1)
        else:
            continue
        
        for idx2, seq in enumerate(df["Peptide"][idx1:]):
            stripped_seq1 = re.sub("\{","", seq)
            stripped_seq1 = re.sub("\}","", stripped_seq1)
            stripped_seq1 = re.sub("\+[0-9]+\.*[0-9]*","", stripped_seq1)
            if stripped_seq1 == stripped_seq:
                modifided_mass = re.findall("[0-9]+\.*[0-9]*", seq)
                modified_mass_sum =  sum(map(float,modifided_mass))
                sequence_dic[stripped_seq].append((query_seq, modified_mass_sum, idx1+idx2))
            
    
    for key in sequence_dic:
        with open(f"./data/{key}.tsv", "w") as f:
            for elements in sequence_dic[key]:
                f.write(f"{elements[0]}\t{elements[1]}\t{elements[2]}\n")
            
    
group_seqeunce()
exit()



def get_SA():
    spec1 = "spectrum path 1"
    spec2 = "spectrum path 2"
    mgf = open("mgf file path")
    mgf_info = mgf.readlines()
    index1 = 0
    index2 = 0   