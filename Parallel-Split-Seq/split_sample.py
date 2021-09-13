import datetime

def get_sample_bc(sample_bc_list):
	sample_bc_dic={}
	with open(sample_bc_list) as input:
		for line in input:
			if line.strip():
				[sample, bc]=line.strip().split()
				sample_bc_dic[bc]=sample
	return sample_bc_dic

def get_prefix(prefix_list):
    # items split by space or tab
    # no space in each item name
    names=[]
    with open(prefix_list) as input:
        for line in input:
            line=line.strip()
            if line:
                names.append(line.split()[0])
    return names

def resplit_bc(output_dir, sample_bc_list, prefix_list):
	### split sample by bc3 ### actually bc1
	# new sample name cannot be the same with any of the original prefix
	wkdir=output_dir+'/parse_bc'
	resplit_dir=wkdir
	fout1={}
	fout2={}

	sample_bc=get_sample_bc(sample_bc_list)
	for bc in sample_bc:
		fout1[bc]=open(resplit_dir+'/'+sample_bc[bc]+'_tag_barcode_R1.fastq','w')
		fout2[bc]=open(resplit_dir+'/'+sample_bc[bc]+'_tag_barcode_R2.fastq','w')

	prefix_list=get_prefix(prefix_list)
	print('prefix_list', prefix_list)
	for prefix in prefix_list:
		print(datetime.datetime.now(), 'Processing prefix...',prefix)
		count={}
		for bc in sample_bc:
			count[bc]=0
		count['OTHERS']=0
		fastq1=wkdir+'/'+prefix+'_tag_barcode_R1.fastq'
		fastq2=wkdir+'/'+prefix+'_tag_barcode_R2.fastq'
		with open(fastq1) as f1, open(fastq2) as f2:
			while True:
				header2 = f2.readline()
				if len(header2)==0:
					break
				bc3=header2.split(':')[0][-6:]

				header1=f1.readline()
				seq1=f1.readline()
				strand1=f1.readline()
				qual1=f1.readline()
				
				seq2=f2.readline()
				strand2=f2.readline()
				qual2=f2.readline()
				
				try:
					fout1[bc3].write(header1+seq1+strand1+qual1)
					fout2[bc3].write(header2+seq2+strand2+qual2)
					count[bc3]+=1
				except:
					count['OTHERS']+=1
		with open(resplit_dir+'/'+prefix+'_resplit_stats.txt','w') as output:
			total_counts=sum(count.values())
			output.write('total_counts\t%d'%total_counts)
			for bc in sample_bc:
				output.write('\n%s\t%.3f'%(sample_bc[bc], count[bc]/total_counts))
			output.write('\nOTHERS\t%.3f'%(count['OTHERS']/total_counts))
            
	for bc in sample_bc:
		fout1[bc].close()
		fout2[bc].close()

	samples=list(sample_bc.values())
	samples=list(set(samples))
	with open(output_dir+'/sample.list','w') as f:
		f.write(samples[0])
		for s in samples[1:]:
			f.write('\n'+s)
            
def resplit_bc(output_dir, sample_bc_list, prefix_list):
        ### split sample by bc3 ### actually bc1
        # new sample name cannot be the same with any of the original prefix
        wkdir=output_dir+'/parse_bc'
        resplit_dir=wkdir
        fout1={}
        fout2={}

        sample_bc=get_sample_bc(sample_bc_list)
        for bc in sample_bc:
                fout1[bc]=open(resplit_dir+'/'+sample_bc[bc]+'_tag_barcode_R1.fastq','w')
                fout2[bc]=open(resplit_dir+'/'+sample_bc[bc]+'_tag_barcode_R2.fastq','w')

        prefix_list=get_prefix(prefix_list)
        print('prefix_list', prefix_list)
        for prefix in prefix_list:
                print(datetime.datetime.now(), 'Processing prefix...',prefix)
                count={}
                for bc in sample_bc:
                        count[bc]=0
                count['OTHERS']=0
                fastq1=wkdir+'/'+prefix+'_tag_barcode_R1.fastq'
                fastq2=wkdir+'/'+prefix+'_tag_barcode_R2.fastq'
                with open(fastq1) as f1, open(fastq2) as f2:
                        while True:
                                header2 = f2.readline()
                                if len(header2)==0:
                                        break
                                #bc3=header2.split(':')[0][-6:]
                                bc=header2.split(':')[0][-6:]

                                header1=f1.readline()
                                seq1=f1.readline()
                                strand1=f1.readline()
                                qual1=f1.readline()
                                
                                seq2=f2.readline()
                                strand2=f2.readline()
                                qual2=f2.readline()
                                
                                try:
                                        fout1[bc].write(header1+seq1+strand1+qual1)
                                        fout2[bc].write(header2+seq2+strand2+qual2)
                                        count[bc]+=1
                                except:
                                        count.setdefault(bc, 0)
                                        count[bc]+=1
                                        #count['OTHERS']+=1
                with open(resplit_dir+'/'+prefix+'_resplit_stats.txt','w') as output:
                        total_counts=sum(count.values())
                        output.write('total_counts\t%d'%total_counts)
                        for bc in sample_bc:
                                output.write('\n%s\t%.3f'%(sample_bc[bc], count[bc]/total_counts))
                        #output.write('\nOTHERS\t%.3f'%(count['OTHERS']/total_counts))
                with open(resplit_dir+'/'+prefix+'_resplit_stats2.txt','w') as output:
                        total_counts=sum(count.values())
                        output.write('total_counts\t%d'%total_counts)
                        for bc in count.keys():
                                output.write('\n%s\t%.3f'%(bc, count[bc]/total_counts))
            
        for bc in sample_bc:
                fout1[bc].close()
                fout2[bc].close()

        samples=list(sample_bc.values())
        samples=list(set(samples))
        with open(output_dir+'/sample.list','w') as f:
                f.write(samples[0])
                for s in samples[1:]:
                        f.write('\n'+s)

            
def split_atac_rna(output_dir, sample, config):
	### split atac and rna ###
	parse_bc_dir=output_dir+'/parse_bc'
	split_dir=output_dir+'/split_sample'
	fastq1=parse_bc_dir+'/'+sample+'_tag_barcode_R1.fastq'
	fastq2=parse_bc_dir+'/'+sample+'_tag_barcode_R2.fastq'
	counts = {'A':0,'R':0,'W':0}


	fout1={'A':open(split_dir+'/'+sample+'_ATAC_barcode_R1.fq','w'), 
		'R':open(split_dir+'/'+sample+'_RNA_barcode_R1.fq','w'),
		'W':open(split_dir+'/'+sample+'_UNK_barcode_R1.fq','w')}
	fout2={'A':open(split_dir+'/'+sample+'_ATAC_barcode_R2.fq','w'), 
		'R':open(split_dir+'/'+sample+'_RNA_barcode_R2.fq','w'),
		'W':open(split_dir+'/'+sample+'_UNK_barcode_R2.fq','w')}
	with open(fastq1) as f1, open(fastq2) as f2:
		while True:
			header2 = f2.readline()
			if len(header2)==0:
				break
			tag=header2.split('_')[0][1]
			header1=f1.readline()
			seq1=f1.readline()
			strand1=f1.readline()
			qual1=f1.readline()	
			seq2=f2.readline()
			strand2=f2.readline()
			qual2=f2.readline()	
			if tag=='A' and len(seq2)<config.atac_mlen2:
				continue
			# remove tag from header
			header1='@'+header1.split('_')[1]
			header2='@'+header2.split('_')[1]
			fout1[tag].write(header1+seq1+strand1+qual1)
			fout2[tag].write(header2+seq2+strand2+qual2)
			counts[tag]+=1
	with open(split_dir+'/'+sample+'_ATAC_pipeline_stats.txt','w') as output:
		output.write('fastq_reads\t'+str(counts['A']))
	with open(split_dir+'/'+sample+'_RNA_pipeline_stats.txt','w') as output:
		output.write('fastq_reads\t'+str(counts['R']))
	with open(split_dir+'/'+sample+'_UNK_pipeline_stats.txt','w') as output:
		output.write('fastq_reads\t'+str(counts['W']))
	with open(split_dir+'/'+sample+'_split_pipeline_stats.txt','w') as output:
		output.write('fastq_reads_total\t'+str(counts['A']+counts['R']+counts['W'])+'\n')
		output.write('fastq_reads_ATAC\t'+str(counts['A'])+'\n')
		output.write('fastq_reads_RNA\t'+str(counts['R'])+'\n')
		output.write('fastq_reads_UNKNOWN\t'+str(counts['W']))

