import pandas as pd
import numpy as np
import re
# import pysptools.distance.SAM
from pysptools.distance import SAM


# modified sequence와 unmodified sequence를 그루핑하는 함수
def make_stripped_seq():
    df = pd.read_csv("chick_modplus_fdr1.tsv", delimiter="\t")
    print(df.shape)
    df = df.loc[:,["SpectrumFile", "Index", "ObservedMW", "Peptide"]]
    print(df.shape)
    
    stripped_seq_list = []
    modified_mass_list = []
    # sequence_dic = {}
    #stripped sequence column을 하나더 만들어 나중에 pd.groupby를 통해 그루핑을 할것
    for idx1, query_seq in enumerate(df["Peptide"]):
        # modifided_mass = re.findall("[0-9]+\.*[0-9]*", query_seq)
        if idx1 % 10000 == 0:
            print(idx1)
        
        stripped_seq = re.sub("\{","", query_seq)
        stripped_seq = re.sub("\}","", stripped_seq)
        stripped_seq = re.sub("(\+|\-)[0-9]+\.*[0-9]*","", stripped_seq)
        stripped_seq_list.append(stripped_seq)

        modifided_mass_plus = re.findall("\+[0-9]+\.*[0-9]+", query_seq)
        modifided_mass_minus = re.findall("\-[0-9]+\.*[0-9]+", query_seq)
        for i in range(len(modifided_mass_plus)):
            modifided_mass_plus[i] = modifided_mass_plus[i][1:]
        for i in range(len(modifided_mass_minus)):
            modifided_mass_minus[i] = modifided_mass_minus[i][1:]

        modified_mass_sum_plus =  sum(map(float,modifided_mass_plus))
        modified_mass_sum_minus =  sum(map(float,modifided_mass_minus))
        modified_mass_sum = modified_mass_sum_plus - modified_mass_sum_minus
        modified_mass_list.append(modified_mass_sum)

    series_stripped_seq = pd.Series(stripped_seq_list, name="seq_stripped")
    series_modified_mass = pd.Series(modified_mass_list, name="modified_mass")

    df = pd.concat([df, series_stripped_seq], axis=1)
    df = pd.concat([df, series_modified_mass], axis=1)
    df.to_csv("with_stripped_seq", sep="\t", index=False)
    
make_stripped_seq()
exit()

def groupby_seq():
    df = pd.read_csv("with_stripped_seq", delimiter="\t")
    df = df.groupby("seq_stripped")
    # for i,a in enumerate(df):
    #     # for b in a:
    #     #     print(b)
    #     print(a)
    #     if i == 10:
    #         break
    print(df.first())
    # print(df.loc[:10])
groupby_seq()
exit()