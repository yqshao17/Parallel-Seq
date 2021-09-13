# samtools merge test_filtered.bam *bam
# input atac bam file for tss enrichment and fragment distribution


import matplotlib as mpl 
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np 
import pysam
import os
import pandas as pd 
from scipy.stats import gaussian_kde
#import seaborn as sns
from collections import Counter

def get_UmiCount(samplebam):
	bamfile = pysam.Samfile(samplebam, "rb")
	UmiCount={}
	for p2_rds in bamfile:
		if p2_rds.is_unmapped: continue
		if p2_rds.mapq<30:continue
		if p2_rds.is_reverse:continue
		else:
			bc=p2_rds.query_name.split(':')[0]
			UmiCount.setdefault(bc,0)
			UmiCount[bc]+=1
	return UmiCount

def get_InsertMat(tssbed, samplebam):	
	InsertMat={}
	def asn_mat(bc,pos,InsertMat,s_int,e_int,ilen,weight):
		InsertMat.setdefault(bc,np.zeros(4000))
		if float(pos)>=s_int and float(pos)<e_int-1 and ilen<1000:
			InsertMat[bc][pos-s_int] += weight
		return InsertMat

	bamfile = pysam.Samfile(samplebam, "rb")
	for line in open(tssbed):
		items=line.strip().split('\t')
		chrx_tss,t1,t2=items[:3]
		tss=int((int(t1)+int(t2))/2)
		s_int=tss-2000
		e_int=tss+2000
		for p2_rds in bamfile.fetch(chrx_tss, max(0,s_int-2000), e_int+2000):
			if p2_rds.is_unmapped: continue
			if p2_rds.mapq<30:continue
			if p2_rds.is_reverse:continue
			else:
				bc=p2_rds.query_name.split(':')[0]
				l_pos = p2_rds.pos+4
				ilen = abs(p2_rds.tlen)-9
				r_pos=p2_rds.pos+abs(p2_rds.tlen)-5
				InsertMat = asn_mat(bc,l_pos,InsertMat,s_int,e_int,ilen,1)
				InsertMat = asn_mat(bc,r_pos,InsertMat,s_int,e_int,ilen,1)
	return InsertMat

def get_TssEnrichmentMat(insert_mat):
	return insert_mat/np.mean(np.hstack((insert_mat[:100],insert_mat[-100:])))

def SmoothMat(mat,k):
        return np.convolve(mat,np.ones(k),'same')/k

def plot_tss_enrich(tssbed, samplebam, prefix):
	UmiCount=get_UmiCount(samplebam)
	InsertMat=get_InsertMat(tssbed,samplebam)	
	TssEnrich=[] # barcode: float
	UMI=[]
	for bc in InsertMat.keys():
		TssEnrich.append(SmoothMat(get_TssEnrichmentMat(InsertMat[bc]),50)[2000])
		UMI.append(UmiCount[bc])

	OverallMat=sum(InsertMat.values())
	OverallTssEnrichMat=get_TssEnrichmentMat(OverallMat)
	OveralSmooth=SmoothMat(OverallTssEnrichMat,50)
	TotalUmi=sum(UMI)

	outdf=pd.DataFrame()
	outdf['barcode']=list(InsertMat.keys())
	outdf['UMI']=UMI
	outdf['TssEnrichScore']=TssEnrich
	outdf.to_csv(prefix+'_TssEnrichUMI.txt', sep='\t', header=True, index=None)

	# plot overalll tss enrichment score
	plt.figure(figsize=(6,4))
	plt.scatter(range(-2000,2000),OverallTssEnrichMat,alpha=0.5,s=4,c='k')
	smooth_mat=np.convolve(OverallMat,np.ones(50),'same')/50
	plt.plot(range(-2000,2000),get_TssEnrichmentMat(smooth_mat),'r')
	plt.xlabel('Distance to TSS center')
	plt.ylabel('TSS enrichment score')
	plt.savefig(prefix+'_Overall_TssEnrichment%.2f_%d.pdf'%(OveralSmooth[2000], TotalUmi))

	# plot TSS enrichment score vs umi counts
	plt.figure(figsize=(6,6))
	x=np.log10(np.array(UMI)+1)
	y=np.array(TssEnrich)
	plt.scatter(x,y,alpha=0.5,s=4)
	plt.xlabel('Log10 UMI counts')
	plt.ylabel('TSS enrichment score')
	plt.title('TSS enrichment score vs UMI counts')
	plt.savefig(prefix+'_TSSEnrichmentScore_UMI.pdf')

	# plot TSS enrichment score with color map
	ix=np.where(~np.isnan(y)&~np.isinf(y))[0]
	x=x[ix]
	y=y[ix]
	plt.figure()
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	plt.scatter(x,y,c=z,s=4,alpha=0.5)
	plt.colorbar()
	y=np.sort(y)
	ylim=y[int(len(y)*0.9)] # show 90 percent cell barcodes
	plt.ylim(0,ylim)
	plt.xlabel("Log10 UMI")
	plt.ylabel("TSS enrichment")
	plt.savefig(prefix+'_TSSEnrichmentScore_colormap.pdf')


def is_sorted_coordinate(header):                                                                                                                                                                                                                                            
	"""
	Check if bam fiel is sorted by read name.
	"""
	if("HD" in header):
		if("SO" in header["HD"]):
			if(header["HD"]["SO"] == "coordinate"):
				return True
	else:
		return False

def TSS_ENRICH(tssbed,samplebam,prefix):
	num_threads=4
	# check sorted
	sortedbam = None
	samfile = pysam.AlignmentFile(samplebam, "rb")
	if not is_sorted_coordinate(samfile.header):
		sortedbam = samplebam[:-4]+'.sorted.bam'
		os.system('samtools sort %s -o %s -@ %d'%(samplebam, sortedbam, num_threads))
	samfile.close()
	
	if sortedbam:
		os.system('samtools index '+sortedbam)
		plot_tss_enrich(tssbed,sortedbam,prefix)
	else:
		os.system('samtools index '+samplebam)
		plot_tss_enrich(tssbed,samplebam,prefix)


def FRAG_DIST(frag_ref,samplebam,prefix):
	bamfile = pysam.Samfile(samplebam, "rb")
	fragL = []
	for p2_rds in bamfile:
		if p2_rds.is_unmapped: continue
		if p2_rds.mapq<30:continue
		if p2_rds.is_reverse:continue
		elif p2_rds.tlen<1000:
			fragL.append(p2_rds.tlen)
	count=Counter(fragL)
	fig=plt.figure(figsize=(5.0,4.0))
	df=pd.read_csv(frag_ref,header=0,index_col=0, sep='\t').iloc[:,1:]
	for i in range(df.shape[1]-1):
		plt.plot(df.index,df.iloc[:,i+1]/1000,color='#DBDBDB',linestyle='-',linewidth=4)
	keys=sorted(list(count.keys()))
	values=[count[k] for k in keys]                                                                                                                                                                                   
	plt.plot(keys,np.array(values)/float(sum(values)),'r')
	plt.xlabel('Read length')
	plt.ylabel('Read counts %') 
	fig.savefig(prefix+'_FragmentDistribution.pdf',bbox_inches='tight')
	plt.close(fig)
	return


def filter_UMI_TssEnrich(bamfile, TssEnrichUMI, cutumi,cuttss,prefix):
	# after caluculating tss enrichment and determine cutoff score
	umi_tss=pd.read_csv(TssEnrichUMI,sep='\t')

	filtered=umi_tss[(umi_tss['UMI']>=cutumi) & (umi_tss['TssEnrichScore']>=cuttss)]
	filtered.index=filtered['barcode']
	d=filtered['UMI'].to_dict()
	output_filename = bamfile[:-3]+ 'filtered.bam'
	inbam=pysam.Samfile(bamfile)
	filtered_bam = pysam.Samfile(output_filename, "wb", template=inbam)
	for read in inbam:
		if read.qname.split(':')[0] in d:
			filtered_bam.write(read)
	inbam.close()
	filtered_bam.close()

	with open(prefix+'_stats.txt','w') as output:
		output.write('total cells:       %d\n'%umi_tss.shape[0])
		output.write('median umi:        %d\n'%np.median(umi_tss['UMI']))
		output.write('median tss score:  %0.2f\n-----\n'%np.median(umi_tss['TssEnrichScore'])) #js_edited_20200219 %d->%0.2f
		output.write('filtered cells:    %d\n'%filtered.shape[0])
		output.write('median umi:        %d\n'%np.median(filtered['UMI']))
		output.write('median tss score:  %0.2f'%np.median(filtered['TssEnrichScore'])) #js_edited_20200219 %d->%0.2f



if __name__ == '__main__':

	import sys
	#samplebam='/Share2/home/zhangqf5/yanqiu/scAR/output/20190823_test/atac_qc/test_filtered.bam'
	#prefix='/Share2/home/zhangqf5/yanqiu/scAR/output/20190823_test/atac_qc/test'
	
	prefix=sys.argv[1]
	
	samplebam=sys.argv[2]
	tssbed=sys.argv[3] 
	#/Share2/home/zhangqf5/yanqiu/scAR/pipelinev2/Data/hg38_mm10_combined_epdnew.bed
	#/Share2/home/zhangqf5/yanqiu/scAR/pipelinev2/Data/human_epdnew_TGfph.bed
	#TSS_ENRICH(tssbed,samplebam,prefix)
	FRAG_DIST('Data/Fragment_length_ratio.txt',samplebam,prefix)

	'''
	umi=sys.argv[2]
	tss_score=sys.argv[3]
	filter_UMI_TssEnrich(bamfile, TssEnrichUMI, cutumi,cuttss,prefix)
	'''
