import pandas as pd
import numpy as np
import re
# import pysptools.distance.SAM
from pysptools.distance import SAM

# string = "UNKNOWN(M1)"
# if "UNKNOWN" in string:
#     print("asd")
# # exit()
def filter_by_annotation(df):
    index_list = []
    for idx, annotation in enumerate(df["ModificationAnnotation"]):
        # print(type(annotation))
        # print(annotation)
        # exit()
        # if annotation is np.nan:
        #     print("asdasd")
        # print(annotation)
        # exit()
        if annotation is np.nan :
            annotation = ""
        
        if "UNKNOWN" in annotation:
            index_list.append(idx)
        elif "LABILE" in annotation:
            index_list.append(idx)
        elif len(annotation.split(" ")) < 2:
            # index_list.append(idx)
            pass
        else:
            index_list.append(idx)
    # print(len(index_list))
    # print(index_list[:100])
    # print(df.shape)
    df.drop(index_list, inplace=True)
    # print(df.shape)
    return df

# modified sequence와 unmodified sequence를 그루핑하는 함수
def make_stripped_seq():
    # df = pd.read_csv("./data/chick_modplus_fdr1.tsv", delimiter="\t")
    # print(df.shape)
    # df = df.loc[:,["SpectrumFile", "Index", "ObservedMW", "Peptide", "ModificationAnnotation","CalculatedMW", "Charge"]]
    # print(df.shape)
    # print(df.columns)
    # df = df[df["Charge"] == 2]
    # print(df.shape)
    # df.to_csv("filtered_charge.tsv", sep="\t", index=False)
    # exit()
    # df = df.loc[(len(df["ModificationAnnotation"]) < 2), :]
    # df = df.drop()
    # df = pd.read_csv("./filtered_charge.tsv", delimiter="\t")
    # df = filter_by_annotation(df)
    # # print(df["ModificationAnnotation"].head(20))
    # # print(df[abs(df["Index"]) < 345].shape)
    # # print(df.shape)
    # # print(df["ModificationAnnotation"])
    # # exit()
    # # print(df.shape)
    # # exit()
    # # print(df["Peptide"].shape)
    # # exit()
    # # df = df.loc[:].
    # df.to_csv("filtered.tsv", sep = "\t" ,index=False)
    # exit()
    df = pd.read_csv("filtered.tsv", delimiter="\t")
    stripped_seq_list = []
    modified_mass_list = []
    print(df.shape)
    # sequence_dic = {}
    #stripped sequence column을 하나더 만들어 나중에 pd.groupby를 통해 그루핑을 할것
    for query_seq in df["Peptide"]:
        # modifided_mass = re.findall("[0-9]+\.*[0-9]*", query_seq)
        # if idx1 % 10000 == 0:
        #     print(idx1)
        # print(df["Peptide"].shape)
        # exit()
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
    print(len(modified_mass_list))
    print(len(stripped_seq_list))
    
    series_stripped_seq = pd.Series(stripped_seq_list, name="seq_stripped")
    series_modified_mass = pd.Series(modified_mass_list, name="modified_mass")
    print(df.shape)
    df = pd.concat([df, series_stripped_seq], axis=1)
    print(df.shape)
    df = pd.concat([df, series_modified_mass], axis=1)
    print(df.shape)
    df.to_csv("with_stripped_seq2", sep="\t", index=False)
    
# make_stripped_seq()

df = pd.read_csv("with_stripped_seq2", delimiter="\t")

series = pd.Series([i for i in df.index], name="idx")
df = df.sort_values(by="seq_stripped")
df = pd.concat([df, series], axis=1)
df.to_csv("sorted.tsv", sep="\t", index=False)
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