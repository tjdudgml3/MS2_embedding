import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
import numpy as np

def draw_histogram():
    cmap = plt.get_cmap('jet')
    low = cmap(0.5)
    medium =cmap(0.25)
    high = cmap(0.7)
    a = cmap(0.9)
    b = cmap(0.05)

    df = pd.read_csv("SA_mean.txt")
    print(f"dataframe shape = {df.shape}")
    df = df.dropna(axis=0)
    print(f"NA sum = {df.isnull().sum()}")
    print(f"dropped dataframe shaep = {df.shape}")
    mod_df = df[df["mod"]==1]
    unmod_df = df[df["mod"]==0]
    unmod_df = unmod_df.sample(n=mod_df.shape[0])
    rand_df = df["rand_SA"].sample(n=mod_df.shape[0])
    # rand_df = []
    print(f"modified dataframe = {mod_df.shape}")
    print(f"unmodified dataframe = {unmod_df.shape}")
    
    plt.hist([unmod_df["SA"],mod_df["SA_mod"], mod_df["SA"], mod_df["SA_same_mod"],rand_df], color=[low, medium, high, a, b], ec='black', bins=15)
    plt.xlabel('SA')
    plt.ylabel('intensity')
    # plt.axis("scaled")
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in [low,medium,high, a, b]]
    labels= ["unmod-unmod","unmod-mod_modSA", "unmod-mod_SA", "mod-mod", "random"]
    plt.legend(handles, labels)
    plt.show()
    
draw_histogram()

def draw_histogram2():
    cmap = plt.get_cmap('jet')
    low = cmap(0.5)
    medium =cmap(0.25)
    high = cmap(0.7)
    a = cmap(0.9)
    b = cmap(0.05)

    df = pd.read_csv("SA_mean.txt")
    print(f"dataframe shape = {df.shape}")
    df = df.dropna(axis=0)
    print(f"NA sum = {df.isnull().sum()}")
    print(f"dropped dataframe shaep = {df.shape}")
    mod_df = df[df["mod"]==1]
    unmod_df = df[df["mod"]==0]
    unmod_df = unmod_df.sample(n=mod_df.shape[0])
    rand_df = df["rand_SA"].sample(n=mod_df.shape[0])
    # rand_df = []
    print(f"modified dataframe = {mod_df.shape}")
    print(f"unmodified dataframe = {unmod_df.shape}")
    # plt.hist([unmod_df["SA"],mod_df["SA_mod"], mod_df["SA"], mod_df["SA_same_mod_with_mod"],rand_df], color=[low, medium, high, a, b], ec='black', bins=15)
    plt.hist([unmod_df["SA"],mod_df["SA_mod"],mod_df["SA_same_mod"]], color=[low, medium, high], ec='black', bins=15)
    plt.xlabel('SA')
    plt.ylabel('intensity')
    # plt.axis("scaled")
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in [low,medium,high]]
    labels= ["unmod-unmod","unmod-mod_modSA", "mod-mod"]
    plt.legend(handles, labels)
    plt.show()

# draw_histogram2()

def draw_scatter():
    df = pd.read_csv("SA_mean.txt")
    df = df[df["mod"]==1]
    # df = df[df["SA"]>0.5]
    df = df.dropna(axis=0)
    # df["SA"] = 2*df["SA"] - 1
    # df["SA_mod"] = 2*df["SA_mod"] - 1
    plt.scatter(df["SA"],df["SA_mod"])
    x = np.linspace(0,1,100)
    plt.xlabel('SA by cosine similarity')
    plt.ylabel('SA by modified cosine similarity')
    plt.plot(x,x, 'r', label='y=x')
    plt.legend()
    plt.show()
    
    
draw_scatter()