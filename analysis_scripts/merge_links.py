import pickle
import pandas as pd
import numpy as np
from glob import glob
import time
import os
import sys
from statsmodels.stats.multitest import fdrcorrection


def merge_links(pattern,gene_name_to_id, gene_id_to_name, gene_tss):
    files = glob(pattern)
    merged_pic = {'peaks':{},
                 'filtered_peaks':{}}
    for f in files:
        pf=pickle.load(open(f,'rb'))
        for gene in pf['peaks']:
            df=pf['peaks'][gene].copy()
            if df.shape[0]==0:
                break
            elif gene in gene_id_to_name:
                gene_id=gene
                gene_name=gene_id_to_name[gene_id]
            elif gene in gene_name_to_id:
                gene_id=gene_name_to_id[gene][0]
                gene_name=gene_id_to_name[gene_id]
            else:
                break
            df['gene_id']=gene_id
            df['gene_name']=gene_name
            df['TSS']=gene_tss[gene_id]
            merged_pic['peaks'][gene_id]=df
            #merged_pic['filtered_peaks'][p]=pf['filtered_peaks'][p]
    '''
    sig_link_df=pd.DataFrame()
    for gene in merged_pic['filtered_peaks']:
        df=merged_pic['filtered_peaks'][gene].copy()
        df['gene_id']=gene
        df['gene_name']=gene_id_to_name[gene]
        df['TSS']=gene_tss[gene]
        df['peak']=list(df.index)
        sig_link_df=sig_link_df.append(df)
    sig_link_df.reset_index(inplace=True)
    '''
    #return sig_link_df
    return merged_pic

def correction_by_all(reg, alpha=0.05):
    reg_fdr={}
    reg_fdr['peaks']={}
    reg_fdr['filtered_peaks']={}
    reg_all=pd.DataFrame()
    index_dict={}
    i=0
    for gene in reg['peaks']:
        df=reg['peaks'][gene]
        df=df[df['zscore']>0]
        index_dict[gene]=range(i, i+df.shape[0])
        reg_all=reg_all.append(df)
        i=i+df.shape[0]
        
    accept, pval_adj=fdrcorrection(reg_all['pval'].values, alpha=alpha)
    reg_all['pval_adj']=pval_adj
    reg_all['accept']=accept
    
    for gene in reg['peaks']:
        df=reg['peaks'][gene]
        df=df[df['zscore']>0]
        df['pval_adj']=reg_all.iloc[index_dict[gene]]['pval_adj'].values
        df['accept']=reg_all.iloc[index_dict[gene]]['accept'].values
        reg_fdr['peaks'][gene]=df
        reg_fdr['filtered_peaks'][gene]=df.loc[df['accept']==True,:]
        
    return reg_fdr, reg_all

def write_filtered_linkes(pkl):
    sig_link_df=pd.DataFrame()
    for gene in pkl['filtered_peaks']:
        df=pkl['filtered_peaks'][gene].copy()
        df['peak']=list(df.index)
        sig_link_df=sig_link_df.append(df)
    sig_link_df.reset_index(inplace=True)
    return sig_link_df

indir=sys.argv[1]
outfile=sys.argv[2]

gene_info_hg38=pickle.load(open('/Share2/home/zhangqf5/yanqiu/library/gene_info_hg38.pkl','rb'))
gene_tss_hg38=pickle.load(open('/Share2/home/zhangqf5/yanqiu/library/gene_tss_hg38.pkl','rb'))
gene_name_to_id=pickle.load(open('/Share2/home/zhangqf5/yanqiu/library/gene_name_to_id_hg.pkl','rb'))
#gene_info_hg19=pickle.load(open('/Share2/home/zhangqf5/yanqiu/library/gene_info_hg19.pkl','rb'))
#gene_tss_hg19=pickle.load(open('/Share2/home/zhangqf5/yanqiu/library/gene_tss_hg19.pkl','rb'))

# merge links
merged_pic = merge_links('%s/*gene_peak_MI_*.pkl'%indir, gene_name_to_id, 
                          gene_info_hg38['gene_id_to_name'], gene_tss_hg38)
pickle.dump(merged_pic, open(outfile+'.tmp','wb'))

# FDR correction
reg_fdr_all, reg_all_df = correction_by_all(merged_pic, alpha=0.05)
pickle.dump(reg_fdr_all,open(outfile, 'wb'))
reg_all_df['peak']=list(reg_all_df.index)
reg_all_df.reset_index(inplace=True)
reg_all_df[reg_all_df['pval']<0.05].to_csv(outfile+'.all.txt', sep='\t', header=True, index=True)
reg_all_df=reg_all_df[reg_all_df['accept']==True]
reg_all_df.to_csv(outfile+'.txt', sep='\t', header=True, index=True)
