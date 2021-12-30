import pandas as pd
import scanpy as sc
import seaborn as sns
import os
import numpy as np
import scipy
import pickle
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus']=False
plt.style.use('ggplot') 
plt.style.use('seaborn-whitegrid')
sc.set_figure_params(dpi=100)
from pairwise_correlation import pairwise_correlation

def normalize(adata, cpm=False, batch_scale=True, target_sum=1e6, chunk_size=20000):
    if cpm:
        sc.pp.normalize_total(adata, target_sum=target_sum)
    if batch_scale:
        adata = batch_scale(adata, chunk_size=chunk_size)
    return adata

def batch_scale(adata, chunk_size=20000):
    """
    Batch-specific scale data
    
    Parameters
    ----------
    adata
        AnnData
    chunk_size
        chunk large data into small chunks
    
    Return
    ------
    AnnData
    """
    for b in adata.obs['batch'].unique():
        idx = np.where(adata.obs['batch']==b)[0]
        scaler = MaxAbsScaler(copy=False).fit(adata.X[idx])
        for i in range(len(idx)//chunk_size+1):
            adata.X[idx[i*chunk_size:(i+1)*chunk_size]] = scaler.transform(adata.X[idx[i*chunk_size:(i+1)*chunk_size]])

    return adata

    
def get_peak(peak):
        if ':' in peak: # chr:start-end
            ch=peak.split(':')[0].split('_')[-1]
            s,e=peak.split(':')[1].split('-')
        elif '_' in peak: # chr_start_end
            ch,s,e=peak.split('_')
        return ch,int(s),int(e)
    
def get_around_peaks(gene_info, rn, an, outprefix, bin_size=250):

    def intersect_features(features, genome_region):
        chosen=[]
        for feature in features:
            ch,s,e=get_peak(feature)
            if ch==genome_region[0]:
                if genome_region[1]<e<genome_region[2] or genome_region[1]<s<genome_region[2]:
                    chosen.append(feature)
        return chosen
    
    def get_gene_id(gene):
        if gene in gene_id_to_name.keys():
            gene_id = gene
        elif gene in gene_name_to_id.keys():
            gene_id = gene_name_to_id[gene][0]
            if len(gene_name_to_id[gene])>1:
                log.write(gene+' has ambiguous gene name. Default choose %s!\n'%gene_id)
                #print(gene+' has ambiguous gene id. Default choose %s!'%gene_id)
        else:
            log.write(gene+' has no valid annotation!\n')
            raise NameError('Not valid gene name: %s!'%gene)
            #print(gene+' has no valid annotation!')
        return gene_id
    
    def around_peaks(gene_id, around_bin, peaks):       
        chrom=gene_id_to_chrom[gene_id].split('_')[-1]
        start=gene_starts[gene_id]  
        end=gene_ends[gene_id]
        strand=gene_strands[gene_id]
        if strand=='-':
            start=end
        start_bin=start-around_bin
        end_bin=start+around_bin
        peaks2=intersect_features(peaks, (chrom, max(0,start_bin), end_bin))
        return peaks2
    
    print(datetime.datetime.now(),'Getting gene around peaks...')
    bin_size=int(bin_size)
    gene_id_to_name=gene_info['gene_id_to_name']
    gene_name_to_id={}
    for gene_id in gene_id_to_name.keys():
        gene_name=gene_id_to_name[gene_id]
        gene_name_to_id.setdefault(gene_name,[])    
        gene_name_to_id[gene_name].append(gene_id)
        gene_name_to_id.setdefault(gene_id.split('.')[0],[])
        gene_name_to_id[gene_id.split('.')[0]].append(gene_id)
    gene_id_to_chrom=gene_info['gene_id_to_chrom']
    gene_starts=gene_info['gene_starts']
    gene_ends=gene_info['gene_ends']
    gene_strands=gene_info['gene_id_to_strand']
    
    log=open(outprefix+'gene_around_peaks_%dk.log'%bin_size,'w')
    gene_around_peaks={}
    for gene in rn.var_names:
            try:
                gene_id=get_gene_id(gene)
            except:
                continue
            gene_around_peaks[gene_id]=[]
            peaks=around_peaks(gene_id, bin_size*1000, an.var_names)
            log.write('%s: %d around peaks\n'%(gene, len(peaks)))
            gene_around_peaks[gene]=peaks

    with open(outprefix+'gene_around_peaks_%dk.pkl'%bin_size, 'wb') as f:
            pickle.dump(gene_around_peaks, f)  
    log.close()
    return gene_around_peaks
    

from scprep.stats import  mutual_information #knnDREMI,
import warnings
warnings.filterwarnings("ignore")


def get_back_peaks(an, peak, correction=True, n_cut=100):
    if not correction:
        back_peaks = np.random.choice(an.var_names, run, replace=True)
    else:
        # choose background peaks with same gc content and accessibility
        #print(1)
        #print(an.var.loc[peak, 'n_cells'])
        #print(an.var['n_cells'][:10].values)
        score1=abs(an.var['n_cells'].values-an.var.loc[peak, 'n_cells'])/an.var.loc[peak, 'n_cells']
        #print(2)
        score2=abs(an.var['gc_content'].values-an.var.loc[peak, 'gc_content'])/an.var.loc[peak, 'gc_content']
        #print(3)
        score=score1*score2
        #print(4)
        back_peaks = an.var_names[np.argsort(score)[:n_cut]]
        #print(5)
    return back_peaks

def get_simi_genes(rna, gene, corr, i, n_cut=100):
    top_ind = np.argsort(corr[i,:])[-n_cut:][::-1]
    genes=rna.var_names[top_ind]
    return genes
    
def get_pval(array, x):
    #array.sort()
    zscore=(x-np.mean(array))/np.std(array)
    pval = scipy.stats.norm.sf(abs(zscore))
    return zscore,pval

def get_gene_peak_dremi(gene, around_peaks, simi_peaks_dict, 
                        Xr,Xp, rna_index_dict, atac_index_dict,
                        correction=True, run=100, pval_cutoff=0.05):
    peaks=pd.DataFrame(index=['MI','zscore','pval'])
    filtered_peaks=pd.DataFrame(index=['MI','zscore','pval'])
    for peak in around_peaks:
            #back_peaks=get_back_peaks(rn, an, gene, peak, correction=correction, run=run)
            back_peaks=simi_peaks_dict[peak][:run]
            back_scores=[mutual_information(get_obs_vector(Xp, p, atac_index_dict),get_obs_vector(Xr, gene, rna_index_dict)) for p in back_peaks]
            dremi=mutual_information(get_obs_vector(Xp, peak, atac_index_dict),get_obs_vector(Xr, gene, rna_index_dict))
            zscore, pval=get_pval(back_scores, dremi)
            peaks[peak]=[dremi, zscore, pval]
            if (pval<pval_cutoff) and (zscore>0):
                filtered_peaks[peak]=[dremi, zscore, pval]
                
    peaks=peaks.T
    filtered_peaks=filtered_peaks.T
    return peaks, filtered_peaks

    
def cal_dremi(rn, an, bin_size, 
              gene_around_peaks,simi_peaks_dict, 
              outprefix, chunk, 
              correction=True, run=100, pval_cutoff=0.05):
        method='MI'
        gene_peak_fit={'peaks':{},'filtered_peaks':{}}
        log=open(outprefix+'gene_peak_%s_%dk_%d.log'%(method,bin_size, chunk),'w')
        
        Xr=rn.X.toarray()
        local_genes=list(rn.var_names)
        rna_index_dict={}
        for i, gene in enumerate(local_genes):
            rna_index_dict[gene]=i

        Xp=an.X.toarray()
        peaks=list(an.var_names)
        atac_index_dict={}
        for i, peak in enumerate(peaks):
            atac_index_dict[peak]=i
        
        del rn
        
        i=0
        for gene in local_genes:
            i+=1
            if i%100==0:
                print(datetime.datetime.now(),'%d genes processed'%i)
            if gene not in gene_around_peaks.keys():
                continue
            around_peaks=list(set(gene_around_peaks[gene])&set(an.var_names))
            if around_peaks:
                peaks, filtered_peaks=get_gene_peak_dremi(gene,around_peaks,simi_peaks_dict, 
                                                          Xr,Xp, rna_index_dict, atac_index_dict,
                                                          correction=correction, run=run, pval_cutoff=pval_cutoff)
                gene_peak_fit['peaks'][gene]=peaks
                gene_peak_fit['filtered_peaks'][gene]=filtered_peaks
                log.write('%s: %d regulatory peaks\n'%(gene, filtered_peaks.shape[0]))
        log.close()
            
        #if outprefix:
        with open(outprefix+'gene_peak_%s_%dk_%d.pkl'%(method,bin_size, chunk), 'wb') as f:
            pickle.dump(gene_peak_fit, f)        
    #return gene_peak_fit

def get_obs_vector(X, feature, index_dict):
    i=index_dict[feature]
    return X[:, i]



            
from glob import glob
def merge_dict(pattern):
    files = glob(pattern)
    merged_pic = {}
    for f in files:
        pf=pickle.load(open(f,'rb'))
        for p in pf:
            merged_pic[p]=pf[p]
    return merged_pic

#'''
import mpi4py.MPI as MPI    
comm = MPI.COMM_WORLD
# the node rank in the whole community
comm_rank = comm.Get_rank()
# the size of the whole community, i.e.,the total number of working nodes in the MPI cluster
comm_size = comm.Get_size()
#'''

if __name__=='__main__':
    import sys
    import datetime
    import random
    rna_file=sys.argv[1]
    atac_file=sys.argv[2]
    inprefix=sys.argv[3]
    outprefix=sys.argv[4] #outdir
    bin_size=int(sys.argv[5])
    batch_norm=False
    
    
    if comm_rank == 0:
        os.makedirs(inprefix, exist_ok=True)
        os.makedirs(outprefix, exist_ok=True)
        #'''
        print(datetime.datetime.now(), 'Read data')
        rna=sc.read_h5ad(rna_file) # rna gene expression
        sc.pp.filter_genes(rna, min_cells=5) #atac.var['n_cells']
        sc.pp.filter_cells(rna, min_genes=100)
        rna.var_names=rna.var['gene_id'].values
        rna.X=rna.layers['logTPM']
        for key in list(rna.layers.keys()):
            del rna.layers[key]  # to reduce memory
            
        atac=sc.read_h5ad(atac_file) # atac matrix with peaks annotation
        atac.X[atac.X>1] = 1
        sc.pp.filter_genes(atac, min_cells=5) #atac.var['n_cells']
        sc.pp.filter_cells(atac, min_genes=100)
        #sc.pp.normalize_total(adata, target_sum=1e6)

        if batch_norm:
            rna = normalize(rna, cpm=False, batch_scale=True)
            atac = normalize(atac, cpm=False, batch_scale=True)
            
        rna.var_names=rna.var['gene_id'].values
        #del rna.uns
        cells=list(set(rna.obs_names)&set(atac.obs_names))
        rna=rna[cells, ]
        atac=atac[cells, ]
        print(datetime.datetime.now(), 'Processed ATAC shape:', atac.shape)
        print(datetime.datetime.now(), 'Processed RNA shape:', rna.shape)
        rna.write_h5ad(outprefix+'rna.h5ad')
        atac.write_h5ad(outprefix+'atac.h5ad')
        
        if os.path.exists('%sgene_around_peaks_%dk.pkl'%(inprefix, bin_size)):
            gene_around_peaks=pickle.load(open('%sgene_around_peaks_%dk.pkl'%(inprefix, bin_size),'rb'))
        else:
            print(datetime.datetime.now(),'Getting around peaks')
            gene_info = pickle.load(open('gene_info.pkl', 'rb'))
            gene_around_peaks = get_around_peaks(gene_info, rna, atac, inprefix, bin_size=bin_size)

        if 'gc_content' not in atac.var.columns:
            print(datetime.datetime.now(),'Calculating GC content of peaks')
            gc=[]
            '''
            import pyfastx
            fa=pyfastx.Fasta('hg38.fa')
            for peak in atac.var_names:
                chrom, start, end=get_peak(peak)
                gc.append(fa[chrom][start:end].gc_content)
            '''    
            from Bio import SeqIO
            from Bio.SeqUtils import GC
            record_dict = SeqIO.to_dict(SeqIO.parse("hg38.fa", "fasta"))
            for peak in atac.var_names:
                chrom, start, end=get_peak(peak)
                gc.append(GC(record_dict[chrom][start:end].seq))
            atac.var['gc_content']=gc
            atac.write_h5ad(outprefix+'atac.h5ad')
            

        if os.path.exists(outprefix+'background_peaks.pkl'):
            simi_peaks_dict=pickle.load(open(outprefix+'background_peaks.pkl', 'rb'))
        else:
            print(datetime.datetime.now(),'Calculating similar peaks')
            simi_peaks_dict={}
            for peak in atac.var_names:
                #print(0)
                peaks=get_back_peaks(atac, peak, correction=True, n_cut=1000)
                #print(1)
                simi_peaks_dict[peak]=peaks
                #print(2)
            pickle.dump(simi_peaks_dict, open(outprefix+'background_peaks.pkl', 'wb'))
    
        all_genes=list(rna.var_names)
        random.shuffle(all_genes)
    
    all_genes = comm.bcast(all_genes if comm_rank == 0 else None, root = 0)
    num_genes = len(all_genes)
    local_offset = np.linspace(0, num_genes, comm_size +1).astype('int')
    local_genes = all_genes[local_offset[comm_rank] :local_offset[comm_rank + 1]]
    
    rna = comm.bcast(rna if comm_rank == 0 else None, root = 0)
    atac = comm.bcast(atac if comm_rank == 0 else None, root = 0)
    gene_around_peaks = comm.bcast(gene_around_peaks if comm_rank == 0 else None, root = 0)
    simi_peaks_dict = comm.bcast(simi_peaks_dict if comm_rank == 0 else None, root = 0)
    
    cal_dremi(rna[:, local_genes], atac, bin_size, 
              gene_around_peaks, simi_peaks_dict, 
              outprefix, comm_rank,
              correction=True, run=100, pval_cutoff=0.05)

    
    
