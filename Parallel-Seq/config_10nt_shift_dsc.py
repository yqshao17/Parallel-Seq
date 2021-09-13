trim_TSO=False
TSO='AAGCAGTGGTATCAACGCAGAGT'
TSO_rev='ACTCTGCGTTGATACCACTGCTT'

a1='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
a2='CTGTCTCTTATACACATCT'

q5=15 # quality filter at 5'
q3=10 # quality filter at 3'
mlen1=61 # min length of read1
mlen2=25 # min length of read2

#method='com'
method='FreeDivergence'

# remove adapter of RNA reads
R1A='TACTGCAGCTGAACCTC' # spacer1 
R2A='AAAAAAAAAAAAAAA'
tag_seq='AGATGTGTATAAGAGACAG'
tag_edit_dist=(2,9)
tag_length=19
bc_edit_dist=(1,1,1)

#umi_bc_len = [11,8,8,6]
#umi_bc_starts = [0,11,36,61]
#tag_start = 67


umi_bc_len = [10,(8,9,10,11),6]
umi_bc_starts = [0,10,(35,36,37,38)]
tag_start=(41,42,43,44)

i5_len = 16

reads_in_cells_thresh=0.92

#mapping
PE=False
ATAC_PE=True
RNA_PE=False
nthreads=4

min_mapq=30
max_flen=1000
atac_mlen2=20

species = 'mm_hg'

if species =='hg':
	tssref='Data/human_epdnew_TGfph.bed'
	refbed='Data/hg38_2k_windows.bed'
elif species =='mm':
	tssref='Data/mouse_epdnew_p44On.bed'
	refbed='Data/mm10_2k_windows.bed'
elif species =='mm_hg':
	tssref='Data/hg38_mm10_combined_epdnew.bed'
	refbed='Data/hg38_mm10_2k_windows.bed'

fragref='Data/Fragment_length_ratio.txt'
cuttss=4
cutumi=100
countby='bin'
