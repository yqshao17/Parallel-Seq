logdir=test/parallel-split/log_RNA
mkdir $logdir
bsub -q Z-ZQF -e $logdir/err -o $logdir/out \
python pipeline.py all_resplit \
--raw_list test/parallel-split/test_ATAC_raw.list \
--sample_bc_list test/parallel-split/test_sample_bc.list \
--output_dir test/parallel-split/ATAC \
--bcset Data/barcodes/bcset0_all \
--keeptype ATAC \
--config config_10nt \
--countby bin \
--genome_dir Data/Ref/INDEX_mm_hg \
--tssref Data/hg38_mm10_combined_epdnew.bed \
--refbed Data/hg38_mm10_2k_windows.bed \
--nthreads 5 \
--logdir $logdir
