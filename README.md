# Installation guide

## Python Dependencies
```
pandas=1.2.4
numpy=1.20.1
scipy=1.6.2
scikit-learn=0.24.1
scanpy=1.8.1
pysam=0.16.0.1
matplotlib=3.3.4
mpi4py=3.0.3
scprep=1.1.0
```
Dependent packages can be installed by ***pip isntall*** conveniently. 

```
git clone https://github.com/yqshao17/Parallel-Seq/
```
Clone the repository to your local directory in a couple of minutes. All scripts can be used directly without installation.

# Parallel-Seq Data Preprocessing
Paired Analysis of RNA and AccessibLe chromatin in the same singLE celL by SEQuencing.

For preprocessing raw sequencing data
```
cd Parallel-Seq
```
## 1. generate genome index
First, you should generate genome index with the required genome fasta and gtf files with one or more species.

```
# example_mkref.sh
python pipeline.py mkref \
    --output_dir  Data/Ref/INDEX_hg_mm \
    --nthreads 16 \
    --genome hg mm \
    --fasta hg38.fa mm10.fa \
    --genes hg38.gtf mm10.gtf
```

Generate bed files as ATAC bin features
```
bedtools makewindows -g Data/hg38.chrom.sizes -w 2000 > Data/hg38_2k_windows.bed
bedtools makewindows -g Data/mm10.chrom.sizes -w 2000 > Data/mm10_2k_windows.bed
```

## 2. prepare input files
1.  prepare fastq files

    For Parallel-Seq, you should convert raw data to fastq files with barcodes in the read name. For Parallel-Split-Seq, skip this step.
```
python dsc_convert.py \
-i test/parallel/raw_data/ATAC_S1_L001 \
-o test/parallel/raw_data/convert_ATAC/test
```

2.  raw_list (sep by tab or space) 
```
# test/parallel/test_RNA_raw.list, test/parallel/test_sample_bc.list
1st column: library name (prefix of output parsed fastq file names, should be unique)
2nd column: library id (8 bytes, last byte is “.”, will be added to cell barcode as batch_id) 
3rd column: R1.fq.gz 
4th column: R2.fq.gz
```

3.  sample_bc_list (sep by tab or space) (optional if barcode1 represent FACS marker)
```
# test/parallel/test_sample_bc.list
1st column: sample name (sample name cannot be the same as any of the original library name)
2nd column: bc1(6nt)
```

## 3. run pipeline

There are different **modes** like pre-process (trimming adapters and parsing barcodes), repslit (split sample by barcode1), post-process (spilt RNA and ATAC, mapping, reads assigmnent and generate count matrix), etc.

Mode: one of "mkref", "resplit", "all_single", "all_batch", "all_resplit", "preproc_single", "preproc_batch", "postproc_single", "postproc_batch". 
```
"all": runs the entire pipeline
"preproc": runs trimming adapter and parsing barcodes
"postproc": runs all steps from spilt RNA and ATAC
"single": input only one sample
"batch": input sample list
"resplit": resplit samples by bc1, input sample_bc_list
```

run pipeline
```
# example_RNA.sh
python pipeline.py all_resplit \
--raw_list test/parallel/test_RNA_raw.list \
--sample_bc_list test/parallel/test_sample_bc.list \
--output_dir test/parallel/RNA \
--genome_dir Data/Ref/INDEX_mm_hg \
--bcset Data/barcodes/bcset0_all \
--config config_10nt \
--keeptype RNA \
--nthreads 5
--logdir test/parallel/log_RNA

# example_ATAC.sh
python pipeline.py all_resplit \
--raw_list test/parallel/test_ATAC_raw.list \
--sample_bc_list test/parallel/test_sample_bc.list \
--output_dir test/parallel/ATAC \
--bcset Data/barcodes/bcset0_all \
--keeptype ATAC \
--config config_10nt \
--countby bin \
--genome_dir Data/Ref/INDEX_mm_hg \
--tssref Data/hg38_mm10_combined_epdnew.bed \
--refbed Data/hg38_mm10_2k_windows.bed \
--nthreads 5 \
--logdir test/parallel/log_ATAC

```

# Analysis tools

Analysis tools for identification of peak-gene associations, mutation-gene associations and eccDNAs in single cells. All scripts should run at analysis_scripts.
```
cd analysis_scripts
```

## 1. Identification of peak-gene associations
An analysis framework based on mutual information of gene expression and chromatin accessibility profiles in the same single-cell to identify peak-gene associations. It is very time and memory consuming for all genes and peaks for a large number of cells. So MPI is used for parallel running.

Input is adata which can be generated by python package Scanpy.

**hg38.fa** is needed for calculating GC content.

**gene_info.pkl** is needed for annotation, which is created when generating genome index in the first step.

```
bsub -q your_queue -n 10 -e err -o out mpirun -np 10 \
    python mpi_mi_regulatory.py \
    rna.h5ad peaks.h5ad \
    regulatory/ regulatory/mi/ 500
```
Then separate link files are merged.
```
python merge_links_new.py regulatory/mi/ regulatory/gene_peak_MI_500k_merged.pkl
```

## 2. Identification of mutations
We identified the major alleles in normal cells as wild-type alleles and the alleles that were different from wild-type one in tumor cells as mutated allele. 
```
bash call_mutation.sh sample1.bam,sample2.bam mutation/ proj_name query_peaks.bed mutation_group.txt
# query_peaks.bed is the region where you want to search mutations, usually the peaks called from ATAC
# mutation_group.txt has two columns with one column as cell barcodes and annother collumn as label of Normal/Tumor

```

## 3. single-cell eccDNA analysis
eccDNAs are identified by [Circle_finder](https://github.com/pk7zuva/Circle_finder/blob/master/circle_finder-pipeline-bwa-mem-samblaster.sh).

```
circle_finder-pipeline-bwa-mem-samblaster.sh 16 hg38.fa R1.fastq R2.fastq 10 outname hg38
```
Then a custom script assigns eccDNAs to single cells.
```
python eccDNA_process.py infile outfile
# infile is *concordant_freq3.2SPLIT-1M.inoneline.txt generated by Circle_finder
```
