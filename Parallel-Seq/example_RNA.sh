logdir=test/parallel/log_RNA
mkdir $logdir
bsub -q Z-ZQF -e $logdir/err -o $logdir/out \
python pipeline.py all_resplit \
--raw_list test/parallel/test_RNA_raw.list \
--sample_bc_list test/parallel/test_sample_bc.list \
--output_dir test/parallel/RNA \
--genome_dir Data/Ref/INDEX_mm_hg \
--bcset Data/barcodes/bcset0_all_dsc \
--config config_10nt_shift_dsc \
--keeptype RNA \
--nthreads 5
