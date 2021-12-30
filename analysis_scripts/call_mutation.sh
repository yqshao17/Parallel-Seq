#indir=$1
sample=$1 #TL210630
outdir=$2 #/Share2/home/zhangqf5/yanqiu/scAR/processed/TLm4
outprefix=$3 #0630
region=$4
meta=$5 #two columns with annotaiton of Normal/Tumor

########### Call mutations ############
# 1. split_cluster_bam
mkdir $outdir/bam
mkdir $outdir/mutation
mkdir $outdir/bam/${outprefix}
mkdir $outdir/mutation/${outprefix}
python split_cluster_bam.py \
    -s $sample \
    -o $outdir/bam/${outprefix} \
    -c $meta \
    --datatype ATAC

# 2.Sort bam file
samtools sort -o $outdir/bam/${outprefix}/Cluster_Tumor_sorted.bam \
    $outdir/bam/${outprefix}/Cluster_Tumor.bam

# 3.samtools mpileup
#$outdir/data_matrix/${outprefix}_peaks/Cluster_Tumor_peaks.narrowPeak
cut -f 1-3  $region \
    > $outdir/mutation/${outprefix}/Tumor_peaks.bed
#/Share/home/zhangqf/usr/bin/
samtools mpileup \
    $outdir/bam/${outprefix}/Cluster_Tumor_sorted.bam \
    -l $outdir/mutation/${outprefix}/Tumor_peaks.bed \
    -f hg38.fa \
    -d 3000 --output-QNAME>$outdir/mutation/${outprefix}/Tumor_peaks_mpileup.txt

# 4.count mutation alleles
python allelecount.py \
    $outdir/mutation/${outprefix}/Tumor_peaks_mpileup.txt \
    $outdir/mutation/${outprefix}/Tumor_peaks_mpileup.pkl \
    allele 10 5


# 5.determine background allele using reads from Normal cells
samtools sort -o $outdir/bam/${outprefix}/Cluster_Normal_sorted.bam \
    $outdir/bam/${outprefix}/Cluster_Normal.bam
samtools mpileup \
    $outdir/bam/${outprefix}/Cluster_Normal_sorted.bam \
    -l $outdir/mutation/${outprefix}/Tumor_peaks_mpileup.pkl.sites \
    -f hg38.fa \
    -d 3000 --output-QNAME>$outdir/mutation/${outprefix}/Normal_mutation_sites_mpileup.txt
python allelecount.py \
    $outdir/mutation/${outprefix}/Normal_mutation_sites_mpileup.txt \
    $outdir/mutation/${outprefix}/Normal_mutation_sites_mpileup_all.pkl \
    all 3 3
