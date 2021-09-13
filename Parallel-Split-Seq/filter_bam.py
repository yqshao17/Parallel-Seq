
# get filtered bam of different species for other analysis
# bsub_script='bsub -q Z-ZQF -n 10 '

from process_atac import split_bam
from multiprocessing import Process
import pysam
import os


def get_filtered_bc(filtered_cell_file):
        bc_list=[]
        with open(filtered_cell_file) as input:
                input.readline()
                for line in input:
                        bc=line.rstrip('\n').split(',')[0]
                        bc_list.append(bc)
        return  bc_list


def get_bc_umi(bc_list, assignment_file):
        bc_umi={}
        for bc in bc_list:
                bc_umi[bc]=[]
        with open(assignment_file) as input:
                input.readline()
                for line in input:
                        items=line.rstrip('\n').split(',')
                        bc=items[0]
                        umi=items[4]
                        if bc in bc_umi.keys():
                                bc_umi[bc].append(umi)
        return bc_umi


def filter_chunk(output_dir,sample,bc_umi,config,chunk=None):
	# filter by read quality and filtered barcode
	# keep the longest read for each umi

	min_mapq=config.min_mapq
	max_flen=config.max_flen

	if chunk is None:
		bamfile = output_dir + '/mapping/'+sample+ '.Aligned.sorted.bam'
		output_filename = output_dir +'/atac_qc/'+sample+ '_filtered.bam'
	else:
		bamfile = output_dir +'/mapping/'+sample+ '.Aligned.sorted.chunk%d.bam' %chunk
		output_filename = output_dir +'/atac_qc/'+sample+ '_filtered.chunk%d.bam' %chunk

	inbam=pysam.Samfile(bamfile)
	filtered_bam = pysam.Samfile(output_filename, "wb", template=inbam)
	
	bc_prev=''
	umi_prev=''
	seqname_prev=''
	for read in inbam:
		rname  = str(read.reference_name)
		if read.is_unmapped or read.mate_is_unmapped or ("chrM" in rname) or read.is_secondary: continue
		if read.is_proper_pair and (abs(read.isize) <= max_flen) and (read.mapq >= min_mapq):
			if read.is_read1:
				seqname=read.qname
				bc=seqname.split('_')[-1].split(':')[0]
				umi=seqname.split(':')[1]
				if bc in bc_umi.keys():
					if umi in bc_umi[bc]:
						if (bc,umi)!=(bc_prev,umi_prev): # keep the first read pair for each umi
							filtered_bam.write(read)
							#filtered_bam.write(inbam.mate(read))
							bc_prev=bc
							umi_prev=umi
							seqname_prev=seqname
			elif read.is_read2: # write read2 if read1 is written
				seqname=read.qname
				if seqname==seqname_prev:
					filtered_bam.write(read)
	filtered_bam.close()
	inbam.close()


def join_bam(output_dir,sample,nthreads):
    filenames = [output_dir +'/atac_qc/'+sample+ '_filtered.chunk%d.bam' %i for i in range(1,nthreads+1)]
    inbam=pysam.Samfile(output_dir + '/mapping/'+sample+ '.Aligned.sorted.bam')
    outbam=pysam.Samfile(output_dir + '/atac_qc/'+sample+ '_filtered.bam', "wb", template=inbam) #,header=inbam.header
    inbam.close()
    for f in filenames:
    	inbam=pysam.Samfile(f)
    	for read in inbam:
    		outbam.write(read)
    	inbam.close()
    outbam.close()
       

def filter_bam(output_dir, sample, config, nthreads):
	""" Gets molecular info for a bam file. Splits the bamfile into 
	nthread chunks and runs in parallel """
	nthreads = int(nthreads)
	split_bam(output_dir, sample, nthreads)

	filtered_cell_file=output_dir+'/DGE_filtered/'+sample+'_cell_metadata.csv'
	assignment_file=output_dir+'/molecule_info/'+sample+'_read_assignment.csv'

	bc_list=get_filtered_bc(filtered_cell_file)
	bc_umi=get_bc_umi(bc_list, assignment_file)

	Pros = []
	for i in range(1,nthreads+1):
		print('Starting thread %d' %i)
		p = Process(target=filter_chunk, args=(output_dir,sample,bc_umi, config, i))
		Pros.append(p)
		p.start()
	for t in Pros:
		t.join()   

	join_bam(output_dir,sample,nthreads)

	for i in range(1,int(nthreads)+1):
		os.remove(output_dir +'/mapping/'+sample+ '.Aligned.sorted.chunk%d.bam' %i)
		os.remove(output_dir +'/atac_qc/'+sample+ '_filtered.chunk%d.bam' %i)



if __name__ == '__main__':
	import sys
	import importlib
	output_dir=sys.argv[1]
	sample=sys.argv[2]
	nthreads=sys.argv[4]
	config = importlib.import_module(sys.argv[3])
	#nthreads=4
	filter_bam(output_dir,sample,config, int(nthreads))

