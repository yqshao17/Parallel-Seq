bsub -q Z-ZQF -e err -o out \
python pipeline.py mkref \
    --output_dir  Data/Ref/INDEX_hg_mm \
    --nthreads 16 \
    --genome hg mm \
    --fasta hg38.fa mm10.fa \
    --genes hg38.gtf mm10.gtf
