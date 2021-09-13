# Parallel-Seq
Paired Analysis of RNA and AccessibLe chromatin in the same singLE celL by SEQuencing.

**Parallel-Seq** and **Parallel-Split-Seq** are two versions of pipeline processing raw sequencing data.

## generate genome index
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

## prepare input files
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

## run pipeline

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

