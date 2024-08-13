import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pickle
import torch
import os

from simulations import simulate_segment, simulate_SNVs, simulate_SNPs

def simulate_data():
    segs = simulate_segment()

    SNV = pd.DataFrame()
    SNP = pd.DataFrame()
    for _, row in segs.iterrows():

        snp = simulate_SNPs(row)
        snv = simulate_SNVs(row)
        
        SNP =  pd.concat([SNP,snp], axis=0)
        SNV =  pd.concat([SNV,snv], axis=0)
        
        SNP = SNP.reset_index(drop=True)
        SNV = SNV.reset_index(drop=True)
        
    return segs, SNV, SNP
    
def plot_data(SNP, SNV, path, save = False):
    sns.set_theme(style="white", font_scale=1.5)
    fig, axes = plt.subplots(2, 3, figsize=(30, 18))

    baf = sns.scatterplot(data=SNP, x="pos", y="baf", s=20, ax=axes[0,0], hue="CN")
    baf.axhline(0.5)
    sns.histplot(data=SNP, x = "baf", ax=axes[0,1], bins = 100, hue="CN",)#multiple="stack"

    dr = sns.scatterplot(data=SNP, x="pos", y="dr", s=20, ax=axes[0,2], hue="CN")
    dr.axhline(1)

    vaf = sns.scatterplot(data=SNV, x="pos", y="vaf", s=20, ax=axes[1,0], hue="CN")
    vaf.axhline(0.5)
    sns.histplot(data=SNV, x = "vaf", ax=axes[1,1], bins = 150, hue="CN",)

    axes[0,0].set_ylim(0,1)
    axes[0,1].set_xlim(0,1)
    axes[1,0].set_ylim(0,1.5)
    axes[1,1].set_xlim(0,1)
    axes[0,2].set_ylim(-1,3) 
    
    if save:
        plt.savefig(path, dpi=300)
    return

def save_data(SNP, SNV, path, name):
    data_input = {}
    data_input['baf'] = torch.tensor(np.expand_dims(SNP.baf,1)).float()
    data_input['dr'] = torch.tensor(np.expand_dims(SNP.dr,1)).float()
    data_input['vaf'] = torch.tensor(np.expand_dims(SNV.vaf,1)).float()
    
    f = open(os.path.join(path,name+'.pkl'),"wb")
    pickle.dump(data_input, f)
    f.close()
    return
    