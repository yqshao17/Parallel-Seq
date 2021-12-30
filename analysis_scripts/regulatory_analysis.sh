rna=$1 # rna adata
atac=$2 # peak adata
outdir=$3
outfile=$4
mpirun -np 20 python mpi_mi_regulatory_new.py $rna $atac $outdir $outdir 500
python merge_links_new.py $outdir $outfile