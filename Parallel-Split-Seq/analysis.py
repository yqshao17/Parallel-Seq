import pandas as pd
import scipy.io as sio
import scipy.interpolate
import numpy as np
import scipy.sparse
import scipy
import gzip
import subprocess
import collections
from collections import defaultdict, Counter
import scipy.sparse as sp_sparse
import warnings
import pickle
import os
#warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import pylab as plt
fsize=12


rc_dict = dict(zip(list('NACGT'),list('NTGCA')))
def reverse_complement(seq):
    return ''.join([rc_dict[s] for s in seq][::-1])

def generate_dge_matrix(df,read_cutoff=10):
    reads_per_cell = df.groupby(df.cell_barcode).size()
    cells = reads_per_cell[reads_per_cell>3]
    all_genes = pd.Series(df.gene.unique()).sort_values()
    all_genes.index = range(len(all_genes))
    gene_dict = dict(zip(all_genes.values,range(len(all_genes.values))))
    cell_dict = dict(zip(cells.index.values,range(len(cells.index.values))))
    rows,cols,vals = [],[],[]
    for bc,g in zip(df.cell_barcode.values,df.gene.values):
        try:
            cell_dict[bc]
        except:
            pass
        else:
            rows.append(cell_dict[bc])
            cols.append(gene_dict[g])
            vals.append(1)
    rows.append(len(cell_dict)-1)
    cols.append(len(gene_dict)-1)
    vals.append(0)
    digital_count_matrix = scipy.sparse.csr_matrix((vals,(rows,cols)),dtype=np.float64)
    thresholded_cells = np.array(digital_count_matrix.sum(1)).flatten()>read_cutoff
    digital_count_matrix = digital_count_matrix[thresholded_cells,:]
    expressed_genes = np.array(digital_count_matrix.sum(0)).flatten()>0
    all_genes = pd.Series(all_genes[expressed_genes])
    digital_count_matrix = digital_count_matrix[:,expressed_genes]
    barcodes = cells.index.values[thresholded_cells]
    return digital_count_matrix,all_genes,barcodes

def barnyard(cell_data,tickstep=10000,s=4,lim=None,ax=None,fig=None):
    species = cell_data.columns[:2]
    colors = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
              (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              'gray']
    #colors = list(sb.color_palette('Set1',n_colors=2)) + ['gray']
    #sb.set_style("white")
    #sb.set_style("ticks")
    
    if ax is None:
        fig = figure(figsize=(3,3))
        ax = fig.add_subplot(111)
    counts1 = cell_data.iloc[:,0]
    counts2 = cell_data.iloc[:,1]
    cell_type1 = counts1>(counts2*9)
    cell_type2 = counts2>(counts1*9)
    mixed_cells = ~(cell_type1|cell_type2)
    plt.scatter(counts1[mixed_cells],
            counts2[mixed_cells],
            color=colors[2],
            s=s,
            label=None)
    plt.scatter(counts1[cell_type2],
            counts2[cell_type2],
            color=colors[0],
            s=s,
            alpha=1,
            label=None)
    plt.scatter(counts1[cell_type1],
            counts2[cell_type1],
            color=colors[1],
            s=s,
            label=None)
    plt.scatter([],[],
            color=colors[0],
            s=10,
            label='%d %s (%0.1f'%(sum(cell_type2),species[1],100*float(sum(cell_type2))/len(cell_type2))+'%)',
           )
    plt.scatter([],[],
            color=colors[1],
            label='%d %s (%0.1f'%(sum(cell_type1),species[0],100*float(sum(cell_type1))/len(cell_type1))+'%)',
            s=10)
    plt.scatter([],[],
            color=colors[2],
            label='%d Mixed (%0.1f'%(sum(mixed_cells),100*float(sum(mixed_cells))/len(mixed_cells))+'%)',
            s=10)

    if lim==None:
        lim = int((counts1+counts2).max()*1.1)
    ax.set_xticks(plt.arange(0,lim,tickstep))
    ax.set_yticks(plt.arange(0,lim,tickstep))
    ax.set_xticklabels(plt.arange(0,lim,tickstep),rotation=90)
    ax.axis([-int(lim/30.),lim,-int(lim/30.),lim])
    ax.set_xlabel('%s UMI Counts' %species[0],fontsize=fsize)
    ax.set_ylabel('%s UMI Counts' %species[1],fontsize=fsize)
    ax.tick_params(labelsize=fsize)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.legend(fontsize=fsize-1,handletextpad=0.025)
    if fig is None:
        return 0
    else:
        return fig,ax

def get_read_threshold(raw_counts, frac_low=0.3, frac_high=0.925, cutoff=100, try_iter=10,window = 20, cellcut=None):
    #raw_counts = raw_counts[raw_counts>10]
    counts = raw_counts.copy()
    counts = counts[counts>cutoff]
    #print(len(raw_counts),len(counts))
    x = np.log10(np.arange(1,len(counts)+1))
    y = np.log10(counts).values
    if cellcut is None:
        f = scipy.interpolate.interp1d(x, y,kind='linear') # smooth the umi-barcode curve
        x_hat = np.linspace(x.min(),x.max(),500)
        y_hat = f(x_hat)
        y_hat = pd.Series(index=x_hat,data=y_hat)
        y_hat_prime = (-y_hat).diff(window).iloc[window:].values # get umi difference between two barocdes in step 20
        threshold = 10**y_hat.iloc[np.argmax(y_hat_prime)]*0.5 # find the largest diff window
        fraction = sum(raw_counts[raw_counts>threshold])/sum(raw_counts)
        #filtered_cells = sum(read_counts>threshold)
        i=1
        while fraction<=frac_low or fraction>frac_high: #filtered_cells>2000:
            if i>try_iter:
                threshold = cutoff
                break
            threshold = 10**y_hat.iloc[y_hat_prime.argsort()[::-1][i]]#*0.5
            fraction = sum(raw_counts[raw_counts>threshold])/sum(raw_counts)
            #print(threshold, fraction)
            #filtered_cells = sum(read_counts>threshold)
            i+=1
    else:
        cellcut=min(cellcut, len(y))
        threshold = 10**(np.sort(y)[-cellcut])
    fraction = sum(raw_counts[raw_counts>threshold])/sum(raw_counts)
    #print(fraction,threshold)
    print(fraction, threshold)
    return threshold

def plot_read_thresh(read_counts,fig=None,ax=None,frac_low=0.3, frac_high=0.925, cutoff=100, try_iter=10, window = 20, cellcut=None):
    read_threshold = get_read_threshold(read_counts,
                                        frac_low=frac_low, 
                                        frac_high=frac_high, 
                                        cutoff=cutoff, 
                                        try_iter=try_iter,
                                        window = window, 
                                        cellcut=cellcut)
    threshold = len(read_counts[read_counts>read_threshold])
    median_umis = read_counts.sort_values(ascending=False)[:threshold].median()
    if ax is None:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
    ax.plot(range(len(read_counts)),
            (read_counts.sort_values(ascending=False)).values,
            color='lightgray',
            linewidth=2)
    ax.plot(range(threshold),
            (read_counts.sort_values(ascending=False)).values[:threshold],
            color='g',
            linewidth=0,marker='.')
    ax.set_xscale('log')
    ax.set_yscale('log')
    _ = ax.set_xlabel('# Barcodes (logscale)')
    _ = ax.set_ylabel('# UMIs (logscale)')
    ax.text(1,10,' n_cells: %d\n read cutoff: %d\n median_umis: %d' %(threshold,read_threshold,median_umis))
    #print(read_threshold)
    if fig is None:
        return read_threshold
    else:
        return fig,ax,read_threshold
    


def generate_single_dge_report(output_dir,genome_dir, datatype, sample_name='',  
                               overwrite=True, analysis_dir=None,
                                frac_low=0.3, frac_high=0.925, 
                               cutoff=100, try_iter=20,
                               window = 20, cellcut=None,
							   rRNA_file='Data/rRNA_hg_mm_dm.fa',
							   mt_rRNA_file='Data/mt_rRNA_hg_mm_dm.fa'):
    
    # Load gene_info dictionary to assign genes to reads
    with open(genome_dir +'/gene_info.pkl', 'rb') as f:
        gene_info = pickle.load(f)

    gene_dict = gene_info['gene_bins']
    exon_gene_start_end_dict = gene_info['genes_to_exons']
    start_dict = gene_info['gene_starts']
    end_dict = gene_info['gene_ends']
    gene_id_to_name = gene_info['gene_id_to_name']
    gene_id_to_genome = gene_info['gene_id_to_genome']

    gene_id_to_name['Intergenic_mm']='Intergenic'
    gene_id_to_name['Ambiguous_mm']='Ambiguous'
    gene_id_to_genome['Intergenic_mm']='mm'
    gene_id_to_genome['Ambiguous_mm']='mm'
    
    gene_id_to_name['Intergenic_hg']='Intergenic'
    gene_id_to_name['Ambiguous_hg']='Ambiguous'
    gene_id_to_genome['Intergenic_hg']='hg'
    gene_id_to_genome['Ambiguous_hg']='hg'
    
    gene_id_to_name['Intergenic_dm']='Intergenic'
    gene_id_to_name['Ambiguous_dm']='Ambiguous'
    gene_id_to_genome['Intergenic_dm']='dm'
    gene_id_to_genome['Ambiguous_dm']='dm'

    # Load the read_assignment file
    df = pd.read_csv(output_dir  + '/molecule_info/' + sample_name + '_read_assignment.csv')
    if datatype=='RNA':
        df=df[(df['gene_name']!='Intergenic')&(df['gene_name']!='Ambiguous')]
        df=df.reset_index()
    #total_reads = df.shape[0]
        
    read_counts = df.groupby('cell_barcode').size().sort_values(ascending=False)
    fig,ax,read_thresh = plot_read_thresh(read_counts,
                                          frac_low=frac_low, 
                                          frac_high=frac_high, 
                                          cutoff=cutoff, 
                                          try_iter=try_iter,
                                          window = window, 
                                          cellcut=cellcut) ###########################
    #print(read_thresh)
    
    digital_count_matrix,all_genes,barcodes = generate_dge_matrix(df,read_cutoff=10)
    
    # gene_df and cell_df
    gene_df = pd.DataFrame()
    gene_df['gene_id'] = all_genes
    gene_df['gene_name'] = all_genes.apply(lambda s:gene_id_to_name[s])
    gene_df['genome'] = all_genes.apply(lambda s:gene_id_to_genome[s])

    species = df.genome.unique()
    species_genes = {}
    species_gene_inds = {}
    species_umi_counts = {}
    species_gene_counts = {}
    for s in species:
        species_genes[s] = all_genes[all_genes.apply(lambda s:gene_id_to_genome[s])==s]
        species_gene_inds[s] = np.where(all_genes.apply(lambda s:gene_id_to_genome[s])==s)[0]
        species_umi_counts[s] = pd.Series(index=barcodes,
                                      data=np.array(digital_count_matrix[:,species_gene_inds[s]].sum(1)).flatten())
        species_gene_counts[s] = pd.Series(index=barcodes,
                                      data=np.array((digital_count_matrix[:,species_gene_inds[s]]>0).sum(1)).flatten())
        
    species_umi_counts = pd.DataFrame(species_umi_counts)
    species_gene_counts = pd.DataFrame(species_gene_counts)
    species_assignments = pd.Series(['multiplet' for i in range(len(barcodes))])
    for s in species:
        species_assignments.loc[np.where((species_umi_counts[s]/species_umi_counts.sum(1))>0.9)] = s

    cell_df = pd.DataFrame()
    cell_df['cell_barcode'] = pd.Series(barcodes)
    cell_df['species'] = species_assignments.values

    cell_df['umi_count'] = np.array(digital_count_matrix.sum(1)).flatten()
    #cell_df['umi_count_50dup'] = cell_df['umi_count'] * 0.5/(1-df.shape[0]/df.counts.sum())
    cell_df['gene_count'] = np.array((digital_count_matrix>0).sum(1)).flatten()
    
    
    # Write unfiltered matrix data
    if not os.path.exists(output_dir + '/DGE_unfiltered/'):
        os.makedirs(output_dir + '/DGE_unfiltered')

    gene_df.to_csv(output_dir + '/DGE_unfiltered/'+  sample_name+ '_genes.csv')
    cell_df.to_csv(output_dir + '/DGE_unfiltered/'+  sample_name+ '_cell_metadata.csv',index=False)
    #sio.mmwrite(output_dir + '/DGE_unfiltered/'+  sample_name+ '_DGE.mtx', digital_count_matrix)
    
    if not os.path.exists(output_dir + '/DGE_unfiltered/'+  sample_name):
        os.makedirs(output_dir + '/DGE_unfiltered/'+  sample_name)
    np.savetxt(output_dir + '/DGE_unfiltered/'+  sample_name+ '/genes.txt', gene_df['gene_id'].values,fmt="%s")
    np.savetxt(output_dir + '/DGE_unfiltered/'+  sample_name+ '/barcodes.txt', cell_df['cell_barcode'].values,fmt="%s")
    sio.mmwrite(output_dir + '/DGE_unfiltered/'+  sample_name+ '/count.mtx', digital_count_matrix.T)
    
    
    # Filter based on automatic cutoff
    valid_cells = np.where(np.array(digital_count_matrix.sum(1)).flatten()>read_thresh)[0]
    digital_count_matrix = digital_count_matrix[valid_cells]
    barcodes = barcodes[valid_cells]
    cell_df = cell_df.iloc[valid_cells]

    #digital_count_matrix.sum()/df.shape[0]
    
    
    # Write filtered matrix data
    if not os.path.exists(output_dir + '/DGE_filtered'):
        os.makedirs(output_dir + '/DGE_filtered')
    gene_df.to_csv(output_dir + '/DGE_filtered/' + sample_name+ '_genes.csv')
    cell_df.to_csv(output_dir + '/DGE_filtered/' + sample_name+ '_cell_metadata.csv',index=False)
    #sio.mmwrite(output_dir + '/DGE_filtered/' + sample_name+ '_DGE.mtx',digital_count_matrix)
    
    if not os.path.exists(output_dir + '/DGE_filtered/'+  sample_name):
        os.makedirs(output_dir + '/DGE_filtered/'+  sample_name)
    np.savetxt(output_dir + '/DGE_filtered/'+  sample_name+ '/genes.txt', gene_df['gene_id'].values,fmt="%s")
    np.savetxt(output_dir + '/DGE_filtered/'+  sample_name+ '/barcodes.txt', cell_df['cell_barcode'].values,fmt="%s")
    sio.mmwrite(output_dir + '/DGE_filtered/'+  sample_name+ '/count.mtx', digital_count_matrix.T)
    

    digital_count_matrix,all_genes,barcodes = generate_dge_matrix(df,read_cutoff=read_thresh)  ########################

    species_genes = {}
    species_gene_inds = {}
    species_umi_counts = {}
    species_gene_counts = {}
    for s in species:
        species_genes[s] = all_genes[all_genes.apply(lambda s:gene_id_to_genome[s])==s]
        species_gene_inds[s] = np.where(all_genes.apply(lambda s:gene_id_to_genome[s])==s)[0]
        species_umi_counts[s] = pd.Series(index=barcodes,
                                      data=np.array(digital_count_matrix[:,species_gene_inds[s]].sum(1)).flatten())
        species_gene_counts[s] = pd.Series(index=barcodes,
                                      data=np.array((digital_count_matrix[:,species_gene_inds[s]]>0).sum(1)).flatten())
        
    species_umi_counts = pd.DataFrame(species_umi_counts)
    species_gene_counts = pd.DataFrame(species_gene_counts)
    species_assignments = pd.Series(['multiplet' for i in range(len(barcodes))])
    for s in species:
        species_assignments.loc[np.where((species_umi_counts[s]/species_umi_counts.sum(1))>0.9)] = s
    species = np.unique(species_assignments.values)
    species = species[species!='multiplet']
    
    sublibraries=None
    if sublibraries is None:
        # Calculate rRNA Percentage:
        kmer_len = 30
        rrna_sense_kmer_dict = {}
        rrna_antisense_kmer_dict = {}
        with open(rRNA_file) as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                if line[0]!='>':
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        rrna_sense_kmer_dict[kmer] = 0
                    line = reverse_complement(line)
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        rrna_antisense_kmer_dict[kmer] = 0
        
        kmer_len = 30
        mt_rrna_sense_kmer_dict = {}
        mt_rrna_antisense_kmer_dict = {}
        with open(mt_rRNA_file) as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                if line[0]!='>':
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        mt_rrna_sense_kmer_dict[kmer] = 0
                    line = reverse_complement(line)
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        mt_rrna_antisense_kmer_dict[kmer] = 0

        def search_kmers(seq,kmer_dict):
            found = False
            for i in range(0,41,10):
                try:
                    kmer_dict[seq[i:i+kmer_len]]
                    found = True
                except:
                    pass
            return found
        
        fastqfile = output_dir + '/split_sample/'+sample_name+'_barcode_R1.fq'
        if os.path.isfile(fastqfile):
            well_counts = {'total_counts':0,
                           'rRNA_sense_counts':0,
                           'rRNA_antisense_counts':0,
                           'mt_rRNA_sense_counts':0,
                           'mt_rRNA_antisense_counts':0}
            #read_lengths = Counter()
            with open(fastqfile) as f:
                for i in range(1000000):
                    header = f.readline()
                    if len(header)==0:
                        break
                    seq = f.readline()[:-1]
                    f.readline()
                    f.readline()
                    
                    well_counts['total_counts'] += 1
                    #read_lengths[len(seq)]+=1
                    if search_kmers(seq,rrna_sense_kmer_dict):
                        well_counts['rRNA_sense_counts'] += 1
                    elif search_kmers(seq,rrna_antisense_kmer_dict):
                        well_counts['rRNA_antisense_counts'] += 1
                    if search_kmers(seq,mt_rrna_sense_kmer_dict):
                        well_counts['mt_rRNA_sense_counts']+= 1
                    elif search_kmers(seq,mt_rrna_antisense_kmer_dict):
                        well_counts['mt_rRNA_antisense_counts'] += 1
        
        #read_len = max(read_lengths.keys())

        rrna_fraction = (well_counts['rRNA_sense_counts']+well_counts['rRNA_antisense_counts'])/well_counts['total_counts']
        mt_rrna_fraction = (well_counts['mt_rRNA_sense_counts']+well_counts['mt_rRNA_antisense_counts'])/well_counts['total_counts']

        
        stat_dict = {}

        with open(output_dir +'/split_sample/'+ sample_name+'_pipeline_stats.txt') as f:
            for line in f:
                k,v = line.strip().split('\t')
                stat_dict[k] = int(v)

        with open(output_dir +'/molecule_info/'+ sample_name+'_pipeline_stats.txt') as f:
            for line in f:
                k,v = line.strip().split('\t')
                stat_dict[k] = int(v)


        stat_dict['Estimated Number of Cells'] = len(barcodes)
        stat_dict['Number of Reads'] = stat_dict['fastq_reads'] #contains all reads, including filtered cells and intergenic reads
        stat_dict['Mean Reads/Cell'] = stat_dict['fastq_reads']/len(barcodes)
        stat_dict['Sequencing Saturation'] = 1-df.shape[0]/df.counts.sum()
        stat_dict['Reads Mapped to rRNA'] = rrna_fraction
        stat_dict['Reads Mapped to mt-rRNA'] = mt_rrna_fraction
        stat_dict['Reads Mapped to genome'] = stat_dict['total_read_count']/stat_dict['fastq_reads']
        stat_dict['Reads Mapped to gene'] = stat_dict['mapped_to_gene']/stat_dict['total_read_count']
        for s in species:
            stat_dict['%s Fraction Reads in Cells' %s] = digital_count_matrix[:,species_gene_inds[s]].sum()/\
                                                         df.query('genome=="%s"' %s, engine='python').shape[0]
            stat_dict['%s Median UMIs/Cell' %s] = np.median(species_umi_counts[s].iloc[np.where(species_assignments==s)])
            #stat_dict['%s Median UMIs/Cell @50%% Dup' %s]  = stat_dict['%s Median UMIs/Cell' %s] * 0.5 /stat_dict['Sequencing Saturation']
            stat_dict['%s Median Genes/Cell' %s] = np.median(species_gene_counts[s].iloc[np.where(species_assignments==s)])
            stat_dict['%s Number of Cells Detected' %s] = sum(species_assignments==s)
            stat_dict['%s Exonic Fraction' %s] = df.loc[np.where(df.cell_barcode.isin(barcodes).values)].query('genome=="%s"' %s).exonic.mean()

        stat_dict['Fraction Reads in Cells'] = digital_count_matrix.sum()/df.shape[0] # umis
    
    stat_catagories = ['Estimated Number of Cells']
    for s in species:
        stat_catagories.append('%s Number of Cells Detected' %s)
    for s in species:
        stat_catagories.append('%s Median UMIs/Cell' %s)    
    #for s in species:
    #    stat_catagories.append('%s Median UMIs/Cell @50%% Dup' %s)
    for s in species:
        stat_catagories.append('%s Median Genes/Cell' %s)
    stat_catagories += ['Mean Reads/Cell',
                        'Number of Reads',
                        'Sequencing Saturation',                        
                        'Reads Mapped to rRNA',
                        'Reads Mapped to mt-rRNA',                       
                        'Reads Mapped to genome',
                        'Reads Mapped to gene',
                        'Fraction Reads in Cells']
    for s in species:
        stat_catagories.append('%s Fraction Reads in Cells' %s)       
    for s in species:
        stat_catagories.append('%s Exonic Fraction' %s)
    
    # Save summary stats to csv
    if overwrite:
        analysis_dir=output_dir+'/analysis'
    else:
        analysis_dir=output_dir+'/'+analysis_dir
        
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    pd.Series(stat_dict).loc[stat_catagories].to_csv(analysis_dir + '/' + sample_name +  '_analysis_summary.csv')
    
    # Subsample reads
    species_read_proportions = df.groupby('genome').size()/df.groupby('genome').size().sum()
    gene_counts_subsampled_df = {}
    umi_counts_subsampled_df = {}
    for s in species:
        seq_depth = species_read_proportions[s] * \
                    stat_dict['Number of Reads']/stat_dict['%s Number of Cells Detected' %s]
        gene_counts_subsampled = {0:0}
        umi_counts_subsampled = {0:0}
        subsample_depths = np.array(list(range(0,int(seq_depth),1000)) + [seq_depth],dtype=int)
        species_df = df.query('genome=="%s"'%s)
        for i in range(1, len(subsample_depths)):
            subsample = subsample_depths[i]
            subsample_fraction = subsample/seq_depth
            sub_sampled_counts = np.random.binomial(species_df.counts.values,subsample_fraction)
            gene_counts_subsampled[subsample] = (species_df[sub_sampled_counts>0]
                                                        .groupby('cell_barcode')
                                                        .gene.apply(lambda x:len(np.unique(x)))
                                                        .reindex(barcodes[np.where(species_assignments==s)])).median()
            umi_counts_subsampled[subsample] = (species_df[sub_sampled_counts>0]
                                                        .groupby('cell_barcode')
                                                        .umi.size()
                                                        .reindex(barcodes[np.where(species_assignments==s)])).median()
        gene_counts_subsampled_df[s] = pd.Series(gene_counts_subsampled).fillna(0)
        umi_counts_subsampled_df[s] = pd.Series(umi_counts_subsampled).fillna(0)
        
    # Generate summary PDF
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0,0.5,0.45,0.5])
    h = 1
    c = 0
    for k in stat_catagories:
        if c < (3*len(species)+3):
            text2write = k+' '*int((34-len(k)))+str(int(np.round(stat_dict[k])))
        else:
            text2write = k+' '*int((34-len(k)))+'%0.3f' %stat_dict[k]
        ax.text(-0.13,h,text2write,fontdict={'family':'monospace'},fontsize=11)
        h-=0.055
        c+=1
    ax.set_axis_off()
    ax = fig.add_axes([0.5,0.65,0.35,0.35])
    _ = plot_read_thresh(read_counts,ax=ax,
                         frac_low=frac_low, 
                         frac_high=frac_high, 
                         cutoff=cutoff, 
                         try_iter=try_iter, 
                         window = window, 
                         cellcut=cellcut)
    ax.set_title(sample_name[:-1])

    if len(species)==2:
        ax = fig.add_axes([1,0.65,0.35,0.35])
        _ = barnyard(species_umi_counts,ax=ax)

    ax = fig.add_axes([0.5,0.1,0.35,0.35])
    for s in species:
        gene_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median Genes per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')

    ax = fig.add_axes([1,0.1,0.35,0.35])
    for s in species:
        umi_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median UMIs per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')
    
    fig.savefig(analysis_dir +'/' + sample_name + '_analysis_summary.pdf',bbox_inches='tight')

    
if __name__=='__main__':
    
    #-c 10 -w 5 -it 20
    import argparse
    import sys
    from split_sample import get_prefix
    from multiprocessing import Process, Manager
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i','--input_dir', help='input dir')
    parser.add_argument('-o', '--output_dir',help='output dir')
    parser.add_argument('-s','--sample', help='sample list', default='sample.list')
    parser.add_argument('-t','--datatype', help='data type: ATAC or RNA')
    parser.add_argument('-fl','--frac_low',type=float, default=0.3)
    parser.add_argument('-fh','--frac_high', type=float, default=0.925)
    parser.add_argument('-c','--cutoff',type=float, help='umi cutoff', default=100)
    parser.add_argument('-it','--try_iter',type=int, help='try knee iterations', default=20)
    parser.add_argument('-w','--window',type=int, help='window', default=20)
    parser.add_argument('-g','--genome_dir', help='genome dir', default='/Share2/home/zhangqf5/yanqiu/scAR_old/split_seq/split-seq-pipeline/INDEX_hg')
    parser.add_argument('--cellcut', help='cells kept', default=None)
    parser.add_argument('--rRNA', help='rRNA.fa', default='Data/rRNA_hg_mm_dm.fa')
    parser.add_argument('--mtrRNA', help='mt_rRNA.fa', default='Data/mt_rRNA_hg_mm_dm.fa')
    args = parser.parse_args()
    
    sample_list=get_prefix(args.input_dir+'/'+args.sample)
    
    '''
    generate_single_dge_report(args.input_dir,
                               args.genome_dir,
                               args.datatype,
                               args.sample+'_'+args.datatype, 
                               overwrite=False, 
                               analysis_dir=args.output_dir,
                               frac_low=args.frac_low, frac_high=args.frac_high, 
                               cutoff=args.cutoff, try_iter=args.try_iter)
    

    '''
    cellcut=args.cellcut
    if cellcut is not None:
        cellcut = int(cellcut)
    manager = Manager()
    Pros = []
    for s in sample_list:
        print('Analyzing sample '+ s)
        p = Process(target=generate_single_dge_report, args=(args.input_dir,
                               args.genome_dir,
                               args.datatype,
                               s+'_'+args.datatype, 
                               False, 
                               args.output_dir,
                               args.frac_low, args.frac_high, 
                               args.cutoff, args.try_iter,
                               args.window,
                               cellcut,
                                args.rRNA, 
                                args.mtrRNA))
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()
        
