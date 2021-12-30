import sys
import re
from collections import Counter
import numpy as np
import pandas as pd
import pickle
import datetime

def parse_mpu(reads, bqs): #, mqs

	# a typical mpileup output line look like (with option: -s -f $refgen)
	# ..+4ACAC.+4ACAC.+2AC => ['.', '.+ACAC', '.+ACAC', '.+AC']
	# .,,-7tttttgtt => ['.', ',', ',-tttttgt', 't']
	# ,$.,,^>. => [',$', '.', ',', ',', '^.']
	# ,.,.*.**.^*. => [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.']
	# ,....,Tn.t => [',', '.', '.', '.', '.', ',', 'T', 'n', '.', 't']
	# A-1N => ['A-N']

	readlist = []
	bqlist = list(bqs)
	#mqlist = list(mqs)

	i = 0		# input pointer in reads
	j = 0		# output pointer in readlist
	while i < len(reads):
		if reads[i] in "ACGTNacgtn.,*#><":
			readlist.append(reads[i])
			i += 1
			j += 1
		elif reads[i] == '$':
			readlist[j-1] += '$'
			i += 1
		elif reads[i] == '^':
			# ^Xa, ^Xa$
			# readlist.append(reads[i:i+3])
			readlist.append(reads[i] + reads[i+2])
			i += 3
			j += 1
		elif reads[i] in '+-':
			# determine length
			digit = re.findall('[\+-](\d+)[ACGTNacgtn*]+',reads[i:])[0]
			readlist[j-1] += reads[i] + reads[i+1+len(digit):i+1+len(digit)+int(digit)]
			i += 1 + len(digit) + int(digit)
		else:
			print('*ERROR* mpileup parser: Unknown char {} in {}[{}]'
				  .format(reads[i], reads, i), file=sys.stderr)
			break

	if len(readlist) != len(bqs):
		print('*ERROR* mpileup parser: length mismatch between BQ string {} '
			  'and reads string {} (breakdown={})'
			  .format(bqs, reads, ':'.join(readlist)), file=sys.stderr)
		return None

	return readlist #, bqlist #, mqlist

def convert_readlist(parse_list, ref):
    merge_dict={".":ref,",":ref,
               "^.":ref,"^,":ref,
               ",$":ref,".$":ref,
               ">":'S',"<":'S',
               "*":'D','#':'D'}
    new_parse_list=[]
    for mut in parse_list:
        if mut in merge_dict:
            new_parse_list.append(merge_dict[mut])
        else:
            new_parse_list.append(mut)
    return new_parse_list

def generate_mut_df(parse_list, cells):
    n_mut=len(set(parse_list))
    n_cells=len(set(cells))
    mut_df=pd.DataFrame(np.zeros((n_cells, n_mut)).astype(int), index=np.unique(cells), columns=np.unique(parse_list))
    for mut, cell in zip(parse_list, cells):
        mut_df.loc[cell, mut]+=1
    return mut_df


if __name__=='__main__':
    import sys
    mpfile = sys.argv[1]
    outpickle = sys.argv[2]
    find = sys.argv[3] # all or allele
    depth_cutoff = int(sys.argv[4]) # 20 filter mutation sites
    cell_cutoff = int(sys.argv[5]) # 10 filter allele

    if find=='all':
        cut_allele = 1
    else:
        cut_allele = 2
    #mpfile='Tumor-CD45-_peaks_mpileup.txt'
    #outpickle='Tumor-CD45-_peaks_mpileup.pkl'
    #outpickle='Tumor-CD45-_peaks_mpileup_new.pkl'
    #mpfile='Normal-CD45-_mutation_sites_mpileup.txt'
    #outpickle='Normal-CD45-_mutation_sites_mpileup_all.pkl'

    #outpickle='Normal-CD45-_peaks_mpileup_all.pkl'
    #outpickle='Normal-CD45-_peaks_mpileup.pkl'

    with open(mpfile) as input:
        kept={}
        n=0
        for line in input:
            #if n>10000:
            #    break
            n+=1
            #if n%1000000==0:
            if n%1000000==0:
                print(datetime.datetime.now(), n, 'processed')
            #if n>1000:
            #    break
            items=line.strip().split('\t')
            chrom, pos, ref, depth, match, qual, names=items
            #mut=['A','G','C','T']
            if int(depth)>=depth_cutoff:
                parse_list = parse_mpu(match.upper().replace(',','.'),qual)
                if parse_list is None:
                    break
                parse_list=convert_readlist(parse_list, ref)
                cells=[el.split(':')[0] for el in names.split(',')]

                # filter mutants pos if there are two mutants with more than 10 reads
                counts = Counter(parse_list)

                '''
                # prefilter
                i=0
                for mut in counts.keys():
                    if mut!='S' and (counts[mut]>=10): # skipped region not considered
                        i+=1
                if i>=2:
                    # filter mutants pos if there are two mutants with more than 10 cells
                    df = generate_mut_df(parse_list, cells)
                    df = df.loc[:,df[df>0].sum(axis=0)>10]
                    #if sum(df[df>0].sum(axis=0)>10)>=2:
                    if df.shape[1]>=2:
                        kept[chrom+':'+pos]=df
                '''

                # not filter mutants for normal samples
                # only filter alleles
                #'''
                i=0
                for mut in counts.keys():
                    if mut!='S' and (counts[mut]>=cell_cutoff): # skipped region not considered, keep allele with more than 10 reads
                        i+=1

                #if i>=1:
                if i>=cut_allele:
                #'''
                #if True:
                    #df = generate_mut_df(parse_list, cells)
                    #df = df.loc[:,df[df>0].sum(axis=0)>10]
                    #new_df = generate_mut_df(parse_list, cells) # slow!!!

                    df = pd.DataFrame({'cells':cells, 'alleles': parse_list})
                    dfg = pd.DataFrame(df.groupby(['cells','alleles']).size())
                    dfg.reset_index(inplace=True)
                    new_df=dfg.groupby(['cells','alleles'])[0].max().unstack()
                    new_df.columns.name = None
                    new_df.index.name = None
                    new_df=new_df.fillna(0)
                    new_df=new_df.astype('int')

                    new_df =new_df.loc[:,new_df[new_df>0].sum(axis=0)>cell_cutoff] # keep allele with more than 10 cells
                    #if new_df.shape[1]>=1:
                    if new_df.shape[1]>=cut_allele:
                        kept[chrom+':'+pos]=new_df

    print(datetime.datetime.now(), 'Writing result')
    pickle.dump(kept, open(outpickle,'wb'))

    if find!='all':
        with open(outpickle+'.sites', 'w') as output:
            for mut in kept.keys():
                chrom=mut.split(':')[0]
                start=mut.split(':')[1]
                output.write(chrom+'\t'+str(int(start)-1)+'\t'+start+'\n')