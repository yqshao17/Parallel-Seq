import pysam
import os
import pandas as pd
import sys
from argparse import ArgumentParser
from glob import glob


'''
split RNA bam files to cluster bam files

usage

python split_cluster_bam.py \
    -s TL201008_new,TL20201006new  \
    -o ../processed/TL_sp_lt/rna_bw2 \
    -c ../processed/TL_sp_lt/sp_lt_seurat_leiden.txt
    
OR

python split_cluster_bam.py -b bam_files \
    -o ../processed/TL_sp_lt/rna_bw2 \
    -c ../processed/TL_sp_lt/sp_lt_seurat_leiden.txt

'''


parser = ArgumentParser(description='call peaks of clusters and get count matrix')
parser.add_argument('-c', help='cluster assignment')
parser.add_argument('-o', help='output directory', default='./')
parser.add_argument('-b', help='input bamfiles', default=None)
#parser.add_argument('-i', help='input directories',default=None)
parser.add_argument('-s', help='samples, sep by comma', default=None)
parser.add_argument('--suffix',action='store_true') # add suffix 0,1,2... for duplicated barcodes from different batches
parser.add_argument('--datatype',default='RNA') 

args = parser.parse_args()

cluster_file=args.c
outdir=args.o
samples = args.s
if samples:
    samples=samples.split(',')
#indir = args.i
bamfiles = args.b

  
os.makedirs(outdir,exist_ok=True)

# read barcode cluster assignment
cluster_assignment=pd.read_csv(cluster_file,sep='\s', header=None, index_col=0)
unique_clusters = cluster_assignment[1].unique()
print('Number of clusters: %d'%len(unique_clusters))
cluster_assignment=cluster_assignment[1].to_dict()



def get_bam_files(s, datatype):
        if datatype=='RNA':
            #pattern=s+'/RNA/mapping/*rmdup.bam'
            pattern=s+'/RNA/mapping/*.Aligned.sorted.bam'
        else:
            pattern=s+'/ATAC/atac_qc/*filtered.sorted.filtered.bam'
        files=glob(pattern)
        return files

    
# generate output bam files
write_bam = {}

if samples:
    files = get_bam_files(samples[0], args.datatype)
elif bamfiles:
    files=bamfiles.split(',')

inbam=pysam.Samfile(files[0])
for c in unique_clusters:
    write_bam[c] = pysam.Samfile(outdir+'/Cluster_%s.bam'%str(c), "wb", template=inbam)
inbam.close()


if samples: # when input directory
    batches = range(len(samples))
    # write output bam files
    for s, b in zip(samples, batches):
        if args.suffix:
            suffix = '-'+str(b)
        else:
            suffix = ''
        files=get_bam_files(s, args.datatype)
        for f in files:
            print(f)
            inbam=pysam.Samfile(f)
            for read in inbam:
                if read.is_unmapped: continue
                #rname  = str(read.reference_name)
                #read.reference_name = rname.split('_')[-1]  # remove hg_ or  mm_
                seqname= read.qname
                bc=seqname.split('_')[-1].split(':')[0]+suffix
                try:
                    c = cluster_assignment[bc]
                    write_bam[c].write(read)
                except:
                    #print(bc, c, '\t'.join(items[:3]+[bc])+'\n')
                    #break
                    pass
            inbam.close()

else: # when input bam files
    batches = range(len(files))
    for f, b in zip(files, batches):
        if args.suffix:
            suffix = '-'+str(b)
        else:
            suffix = ''
    
        print(f)
        inbam=pysam.Samfile(f)
        for read in inbam:
                if read.is_unmapped: continue
                #rname  = str(read.reference_name)
                #read.reference_name = rname.split('_')[-1]  # remove hg_ or  mm_
                seqname= read.qname
                bc=seqname.split('_')[-1].split(':')[0]+suffix
                try:
                    c = cluster_assignment[bc]
                    write_bam[c].write(read)
                except:
                    #print(bc, c, '\t'.join(items[:3]+[bc])+'\n')
                    #break
                    pass
        inbam.close()


for c in unique_clusters:
    write_bam[c].close()
