# first filter samples bam by filtered barcodes, keep unique umi reads
# merge bam files
# calculate fragment distribution and tss enrichment
# filter bam by umi and tss enrichment score

from FragDist_TSSEnrich import TSS_ENRICH,FRAG_DIST,filter_UMI_TssEnrich
from filter_bam import filter_bam
from split_sample import get_prefix
from peak_count_matrix import summit_extend, get_count_matrix

from argparse import ArgumentParser
import os 

parser = ArgumentParser(description='atac bam filtering')
parser.add_argument('mode', help='filter_bc or plot_umi_tss or  filter_umi_tss')

parser.add_argument('-o', '--outdir', help='output directory of filtered bam file', default='./')
parser.add_argument('-t', '--threads', help='number of threads', default='4')

parser.add_argument('-i', '--inbam', help='bam file for umi_tss filtering', default='input bam file')
parser.add_argument('-r', '--tssref', help='tss reference')
parser.add_argument('-n', '--name', help='output name prefix')
parser.add_argument('-f', '--fragref', help='fragment length ratio reference')

parser.add_argument('-cr', '--cutref', help='info of UMI and TSS enrichment score')
parser.add_argument('-ct', '--cuttss', help='cutoff of tss score')
parser.add_argument('-cu', '--cutumi', help='cutoff of UMI')
parser.add_argument('--countby', help='count by peak or bin or promoter', default='peak')
parser.add_argument('--refbed', help='required if countby bin or promoter', default=None) 
# Data/hg38_10k_windows.bed,Data/hg38_Gencode_promoter.bed


args = parser.parse_args()

mode=args.mode.lower()
nthreads=int(args.threads)
outdir=args.outdir

bsub_script='bsub -q Z-ZQF -e atac.err -o atac.out -n %d' % nthreads

if mode=='filter_bc':
    filter_bam(outdir,args.name,config,nthreads)

elif mode=='plot_umi_tss':
    # -r --tssref
    # -f --fragref
    # -i --inbam
    # -o, -n
    FRAG_DIST(args.fragref,args.inbam,outdir+'/'+args.name)
    TSS_ENRICH(args.tssref,args.inbam,outdir+'/'+args.name)

elif mode=='filter_umi_tss':
    # -i, -o, -n
    # -cr, -ct, -cu
    filter_UMI_TssEnrich(args.inbam, args.cutref, int(args.cutumi), float(args.cuttss), outdir+'/'+args.name)

    
elif mode=='count_matrix':
    # -o, -n, i
    bedfile=outdir+'/'+args.name+'_shift.bed'
    #if (not os.path.exists(bedfile)) or (os.path.exists(bedfile) and os.path.getsize(bedfile)==0):
    os.system("python snap_prebed.py -i %s -o %s -m 30 -t 4 -f 20000 -e 0"%(args.inbam,bedfile))
    if args.countby=='bin' or args.countby=='promoter':
        get_count_matrix(args.refbed, 
                         "%s/%s_shift.bed"%(outdir, args.name), 
                         outdir, args.name+'_'+args.countby)
    elif args.countby=='peak':
        os.system("""bash bash/macs2_py2.sh {0}/{1}_shift.bed {0}/{1}""".format(outdir, args.name))
        summit_extend("%s/%s_summits.bed"%(outdir, args.name), "%s/%s_summits_extend.bed"%(outdir, args.name), 500)
        get_count_matrix("%s/%s_summits_extend.bed"%(outdir, args.name), 
                         "%s/%s_shift.bed"%(outdir, args.name), 
                         outdir, args.name+'_'+args.countby)

