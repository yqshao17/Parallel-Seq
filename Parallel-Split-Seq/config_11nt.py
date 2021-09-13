trim_TSO=False
TSO='AAGCAGTGGTATCAACGCAGAGT'
TSO_rev='ACTCTGCGTTGATACCACTGCTT'
a1="CTGTCTCTTATACACATCT"
a2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
q5=15 # quality filter at 5'
q3=10 # quality filter at 3'
mlen1=25 # min length of read1
mlen2=86 # min length of read2

method='FreeDivergence'

R1A='TACTGCAGCTGAACCTC' # spacer1 
tag_seq='AGATGTGTATAAGAGACAG'
tag_edit_dist=(2,9)
tag_length=19
bc_edit_dist=(1,1,1)

umi_bc_len = [11,8,8,6]
umi_bc_starts = [0,11,36,61]
tag_start = 67

#umi_bc_len = [10,8,8,6]
#umi_bc_starts = [0,10,35,60]
#tag_start=66

reads_in_cells_thresh=0.92

#mapping
PE=False
ATAC_PE=True
RNA_PE=False
nthreads=4

min_mapq=30
max_flen=1000
atac_mlen2=20

cuttss=6
cutumi=100

