# split standard 10x atac fastq files by I1 to different samples

import sys
import gzip
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--in_prefix','-i', help='input directory')
parser.add_argument('--outdir','-o', help='output directory')

args = parser.parse_args()
output_dir = args.outdir
in_prefix = args.in_prefix

os.makedirs(output_dir,  exist_ok = True)


R1=in_prefix+'_R1_001.fastq.gz' # Read 1
R2=in_prefix+'_R2_001.fastq.gz' # i5 index as barcode3
R3=in_prefix+'_R3_001.fastq.gz' # Read 2
I1=in_prefix+'_I1_001.fastq.gz' # i7 index as sample id
fq1 = gzip.open(R1,'rb')
fq2 = gzip.open(R2,'rb')
fq3 = gzip.open(R3,'rb')
fqi = gzip.open(I1,'rb')

files1 = open(output_dir + '/R1_001.fastq' ,'w')
files2 = open(output_dir + '/R2_001.fastq' ,'w')
    
n=0
while True:
#while n<1000000:
    if n%1000000==0:
        print('Processed ',n)
    seqnamei = fqi.readline().decode("utf-8")
    if not seqnamei:
        break
    seqi = fqi.readline().decode("utf-8")
    strandi = fqi.readline().decode("utf-8")
    quali = fqi.readline().decode("utf-8")
    
    seqname1 = fq1.readline().decode("utf-8")
    seq1 = fq1.readline().decode("utf-8")
    strand1 = fq1.readline().decode("utf-8")
    qual1 = fq1.readline().decode("utf-8")
    
    seqname2 = fq2.readline().decode("utf-8")
    seq2 = fq2.readline().decode("utf-8")
    strand2 = fq2.readline().decode("utf-8")
    qual2 = fq2.readline().decode("utf-8")
    
    seqname3 = fq3.readline().decode("utf-8")
    seq3 = fq3.readline().decode("utf-8")
    strand3 = fq3.readline().decode("utf-8")
    qual3 = fq3.readline().decode("utf-8")

    i7 = seqi.strip()
    i5 = seq2.strip()
    
    seqname1=seqname1.strip()+'+'+i5+'\n'
    seqname3=seqname3.strip()+'+'+i5+'\n'
    
    files1.write(seqname1+seq1+strand1+qual1)
    files2.write(seqname3+seq3+strand3+qual3)
    n+=1
    
    
files1.close()
files2.close()
fq1.close()
fq2.close()
fq3.close()
fqi.close()

with open(output_dir+'/stats.txt', 'w') as output:
    output.write('Total counts: '+str(n))

os.system('bsub -q Z-ZQF gzip %s/R1_001.fastq' %output_dir)
os.system('bsub -q Z-ZQF gzip %s/R2_001.fastq' %output_dir)
