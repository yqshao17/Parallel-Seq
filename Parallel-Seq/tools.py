# miscellaneous tools
import os
import subprocess
import sys

import pandas as pd
from collections import defaultdict
import gzip
from numpy import unique
import numpy as np
import pickle
import datetime
import pysam
#import HTSeq


#PATH = './'
PATH = os.path.dirname(__file__)
HOME = os.path.expanduser('~')

STAR_PATH = os.path.join(HOME, 'split_seq_reqs', 'bin', 'STAR')
if not os.path.exists(STAR_PATH):
    STAR_PATH = 'STAR'

SAMTOOLS_PATH = os.path.join(HOME, 'split_seq_reqs', 'bin', 'samtools')
if not os.path.exists(SAMTOOLS_PATH):
    SAMTOOLS_PATH = 'samtools'

CUTADAPT_PATH = 'cutadapt'

def get_rawlist(filename):
    # items split by space or tab
    # no space in each item name
    names=[]
    with open(filename) as input:
        for line in input:
            line=line.strip()
            if line:
                names.append(line.split())
    return names
    
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

def download_genome(genome_dir, ref='hg19'):
    """
    Downloads the hg19 reference genome...
    """
    # TODO: find the hg19 genome???
    
def make_combined_genome(species, fasta_filenames, output_dir):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create a combined fasta file with species names added to the start of each chromosome name
    cur_fa = fasta_filenames[0]
    cur_species = species[0]
    if fasta_filenames[0].split('.')[-1]=='gz':
        command = """gunzip -cd {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' > {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
    else:
        command = """cat {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' > {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
    rc = subprocess.call(command, shell=True)
    
    for i in range(1,len(species)):
        cur_fa = fasta_filenames[i]
        cur_species = species[i]
        if fasta_filenames[0].split('.')[-1]=='gz':
            command = """gunzip -cd {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' >> {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
        else:
            command = """cat {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' >> {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
        rc = subprocess.call(command, shell=True)
        
def split_attributes(s):
    """ Returns a dictionary from string of attributes in a GTF/GFF file
    """
    att_list = s[:-1].split('; ')
    att_keys = [a.split(' ')[0] for a in att_list]
    att_values = [' '.join(a.split(' ')[1:]) for a in att_list]
    return dict(zip(att_keys,att_values))

def get_attribute(s,att):
    att_value = ''
    try:
        att_value = split_attributes(s)[att].strip('"')
    except:
        att_value = ''
    return att_value

def parse_gtf(gtf, attributes):
    # get columns or attributes from gtf file
    names = ['Chromosome','Source','Feature','Start','End','Score','Strand','Frame','Attributes']
    if type(gtf)==str:
        gtf=pd.read_csv(gtf,sep='\t',header=None,index_col=None, comment='#')
        gtf.columns=names
    gtf_attr=pd.DataFrame()
    for attr in attributes:
        if attr in names:
            gtf_attr[attr]=gtf[attr]
        else:
            try:
                gtf_attr[attr]=gtf['Attributes'].apply(lambda s: get_attribute(s,attr))
            except:
                print('ValueError: Attribute Name not valid!')
    return gtf_attr

def make_gtf_annotations(species, gtf_filenames, output_dir, splicing):
    splicing = splicing=='True'

    # Load the GTFs
    names = ['Chromosome',
         'Source',
         'Feature',
         'Start',
         'End',
         'Score',
         'Strand',
         'Frame',
         'Attributes']

    gtfs = {}
    for i in range(len(species)):
        s = species[i]
        filename = gtf_filenames[i]
        gtfs[s] = pd.read_csv(filename,sep='\t',names=names,comment='#',engine='python')
    
    # TODO: allow users to specify the gene biotypes that they want to keep
    # For now we keep the following
    gene_biotypes_to_keep = ['protein_coding',
                             'lincRNA',
                             'antisense',
                             'IG_C_gene',
                             'IG_C_pseudogene',
                             'IG_D_gene',
                             'IG_J_gene',
                             'IG_J_pseudogene',
                             'IG_V_gene',
                             'IG_V_pseudogene',
                             'TR_C_gene',
                             'TR_D_gene',
                             'TR_J_gene',
                             'TR_J_pseudogene',
                             'TR_V_gene',
                             'TR_V_pseudogene']
    if splicing:
        # Generate a combined GTF with only the gene annotations
        gtf_gene_combined = gtfs[species[0]].query('Feature=="gene"')
        gtf_gene_combined.loc[:,'Chromosome'] = species[0] + '_' + gtf_gene_combined.Chromosome.apply(lambda s:str(s))
        for i in range(1,len(species)):
            gtf_gene_combined_temp = gtfs[species[i]].query('Feature=="gene"')
            gtf_gene_combined_temp.loc[:,'Chromosome'] = species[i] + '_' + gtf_gene_combined_temp.Chromosome.apply(lambda s:str(s))
            gtf_gene_combined = pd.concat([gtf_gene_combined,gtf_gene_combined_temp])
        gene_biotypes = gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_biotype'))
        #gtf_gene_combined = gtf_gene_combined.iloc[np.where(gene_biotypes.isin(gene_biotypes_to_keep).values)]
        gtf_gene_combined.index = range(len(gtf_gene_combined))
        gtf_gene_combined.to_csv(output_dir + '/genes.gtf',sep='\t',index=False)
    
    # Generate a combined GTF with only the exon annotations
    gtf_exon_combined = gtfs[species[0]].query('Feature=="exon"')
    gtf_exon_combined.loc[:,'Chromosome'] = species[0] + '_' + gtf_exon_combined.Chromosome.apply(lambda s:str(s))
    for i in range(1,len(species)):
        gtf_exon_combined_temp = gtfs[species[i]].query('Feature=="exon"')
        gtf_exon_combined_temp.loc[:,'Chromosome'] = species[i] + '_' + gtf_exon_combined_temp.Chromosome.apply(lambda s:str(s))
        gtf_exon_combined = pd.concat([gtf_exon_combined,gtf_exon_combined_temp])
    gene_biotypes = gtf_exon_combined.Attributes.apply(lambda s: get_attribute(s,'gene_biotype'))
    #gtf_exon_combined = gtf_exon_combined.iloc[np.where(gene_biotypes.isin(gene_biotypes_to_keep).values)]
    gtf_exon_combined.index = range(len(gtf_exon_combined))
    gtf_exon_combined.to_csv(output_dir + '/exons.gtf',sep='\t',index=False)
    
    if not splicing:
        gtf_gene_combined = gtf_exon_combined.copy(deep=True)
        gtf_gene_combined['Feature'] = 'gene'
        gtf_gene_combined.to_csv(output_dir + '/genes.gtf',sep='\t',index=False)
    # Get locations of genes. We are using the longest possible span of different transcripts here
    gtf_gene_combined.loc[:,'gene_id'] = gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id'))
    gene_starts = gtf_gene_combined.groupby('gene_id').Start.apply(min)
    gene_ends = gtf_gene_combined.groupby('gene_id').End.apply(max)
    chroms = gtf_gene_combined.groupby('gene_id').Chromosome.apply(lambda s:list(s)[0])
    strands = gtf_gene_combined.groupby('gene_id').Strand.apply(lambda s:list(s)[0])
    
    gtf_dict_stepsize = 10000
    # Create a dictionary for each "bin" of the genome, that maps to a list of genes within or overlapping
    # that bin. The bin size is determined by gtf_dict_stepsize.
    starts_rounded = gene_starts.apply(lambda s:np.floor(s/gtf_dict_stepsize)*gtf_dict_stepsize).values
    ends_rounded = gene_ends.apply(lambda s:np.ceil(s/gtf_dict_stepsize)*gtf_dict_stepsize).values
    gene_ids = gene_starts.index
    start_dict = gene_starts.to_dict()
    end_dict = gene_ends.to_dict()
    gene_dict = defaultdict(list)
    for i in range(len(gene_starts)):
        cur_chrom = chroms[i]
        cur_strand = strands[i]
        cur_start = int(starts_rounded[i])
        cur_end = int(ends_rounded[i])
        cur_gene_id = gene_ids[i]
        for coord in range(cur_start,cur_end+1,gtf_dict_stepsize):
            if not (cur_gene_id in gene_dict[cur_chrom + ':' +  str(coord)]):
                gene_dict[cur_chrom + ':' +  str(coord)+':'+cur_strand].append(cur_gene_id)
                
    # Create a dictionary from genes to exons
    exon_gene_ids = gtf_exon_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')).values
    exon_starts = gtf_exon_combined.Start.values
    exon_ends = gtf_exon_combined.End.values
    
    exon_gene_start_end_dict = defaultdict(dict)
    for i in range(len(exon_gene_ids)):
        cur_gene_id = exon_gene_ids[i]
        cur_exon_start = exon_starts[i]
        cur_exon_ends = exon_ends[i]
        exon_gene_start_end_dict[cur_gene_id][cur_exon_start] = cur_exon_ends
        
    gene_id_to_gene_names = dict(zip(gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')),
                                     gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_name'))))
    gene_id_to_genome = dict(zip(gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')),
                                 gtf_gene_combined.Chromosome.apply(lambda s:s.split('_')[0])))
    gene_id_to_strand = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                 gtf_gene_combined.Strand.values))
    gene_id_to_chrom = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                 gtf_gene_combined.Chromosome.values))
    gene_id_to_biotype = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                  gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_biotype')).values))

    # add intergenic and ambiguous region for ATAC statistics
    for s in species:
        gene_id_to_gene_names['Intergenic_'+s]='Intergenic'
        gene_id_to_gene_names['Ambiguous_'+s]='Ambiguous'
        gene_id_to_genome['Intergenic_'+s]=s
        gene_id_to_genome['Intergenic_'+s]=s

    #Save dictionary with gene info
    gene_info = {'gene_bins':gene_dict,
                 'genes_to_exons':exon_gene_start_end_dict,
                 'gene_starts': start_dict,
                 'gene_ends': end_dict,
                 'gene_id_to_name': gene_id_to_gene_names,
                 'gene_id_to_genome':gene_id_to_genome,
                 'gene_id_to_chrom':gene_id_to_chrom,
                 'gene_id_to_strand':gene_id_to_strand,
                 'gene_id_to_biotype':gene_id_to_biotype
                }
    
    with open(output_dir+ '/gene_info.pkl', 'wb') as f:
        pickle.dump(gene_info, f, pickle.HIGHEST_PROTOCOL)
        
def generate_STAR_index(output_dir, nthreads,genomeSAindexNbases,splicing):
    splicing = (splicing=='True')
    if splicing:
        star_command = """STAR  --runMode genomeGenerate --genomeDir {0} --genomeFastaFiles {0}/genome.fa --sjdbGTFfile {0}/exons.gtf --runThreadN {1} --limitGenomeGenerateRAM 24000000000 --genomeSAindexNbases {2}""".format(output_dir, nthreads, genomeSAindexNbases)
    else:
        star_command = """STAR  --runMode genomeGenerate --genomeDir {0} --genomeFastaFiles {0}/genome.fa --runThreadN {1} --limitGenomeGenerateRAM 24000000000 --genomeSAindexNbases {2}""".format(output_dir, nthreads, genomeSAindexNbases)
    rc = subprocess.call(star_command, shell=True)
    return rc

bases = list('ACGT')
def convert_degen_seq_to_list(seq):
    """Uses recursion to convert a degenerate sequence to a list
    For example: AGGN -> [AGGA, AGGC, AGGG, AGGT]"""

    seq_list = []
    N_pos = seq.find('N')
    if N_pos>=0:
        for b in bases:
            seq_list += convert_degen_seq_to_list(seq[:N_pos] + b + seq[N_pos+1:])
    else:
        seq_list.append(seq)
    return seq_list

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def get_min_edit_dists(bc,edit_dict,max_d=3):
    """Returns a list of nearest edit dist seqs
    Input 8nt barcode, edit_dist_dictionary
    Output <list of nearest edit distance seqs>, <edit dist>"""
    bc_matches = edit_dict[0][bc]
    edit_dist = 0
    if (len(bc_matches)==0) and (max_d>=1):
        edit_dist+=1
        bc_matches = edit_dict[1][bc]
    if (len(bc_matches)==0) and (max_d>=2):
        edit_dist+=1
        bc_matches = edit_dict[2][bc]
    if (len(bc_matches)==0) and (max_d>=3):
        edit_dist+=1
        bc_matches = edit_dict[3][bc]
    return bc_matches,edit_dist

def trim(output_dir, sample, fq1, fq2, config, nthreads):
    """cut adapter and quality filtering
    """
    trimmed_dir=output_dir+'/trimmed'
    nthreads=int(nthreads)
    q5=config.q5
    q3=config.q3
    mlen1=config.mlen1
    mlen2=config.mlen2
    a1=config.a1
    a2=config.a2 
    if config.trim_TSO==False:
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1},{2} -m {3}:{4} -a {5} -A {6} -o {7}/{8}_R1.fq.gz -p {7}/{8}_R2.fq.gz {9} {10} 1>{7}/{8}.log""".format(nthreads,q5,q3,mlen1,mlen2,a1,a2,trimmed_dir,sample,fq1,fq2), shell=True)
    if config.trim_TSO==True:
        a1=config.TSO
        a2=config.TSO_rev
        subprocess.call(CUTADAPT_PATH + """ --pair-filter=first --discard-untrimmed -j {0} -q {1},{2} -g {3} -A {4} -o {5}/{6}_tmp_R1.fq.gz -p {5}/{6}_tmp_R2.fq.gz {7} {8} 1>{5}/{6}_tmp.log""".format(nthreads,q5,q3,a1,a2,trimmed_dir,sample,fq1,fq2), shell=True)
        
        a1=config.a1
        a2=config.a2
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1},{2} -m {3}:{4} -a {5} -A {6} -o {7}/{8}_R1.fq.gz -p {7}/{8}_R2.fq.gz {7}/{8}_tmp_R1.fq.gz {7}/{8}_tmp_R2.fq.gz 1>{7}/{8}.log""".format(nthreads,q5,q3,mlen1,mlen2,a1,a2,trimmed_dir,sample,fq1,fq2), shell=True)       
    return

def trim_RNA(output_dir, sample, config, nthreads):
    trimmed_dir=output_dir+'/split_sample'
    nthreads=int(nthreads)
    q5=config.q5
    q3=config.q3
    mlen1=config.mlen1
    mlen2=config.mlen2
    a1=config.R1A
    a2=config.R2A
    fq1c=trimmed_dir+'/'+sample+'_barcode_R1.fq'
    fq2c=trimmed_dir+'/'+sample+'_barcode_R2.fq'
    fq1=trimmed_dir+'/'+sample+'_barcode_noretrim_0_R1.fq'
    fq2=trimmed_dir+'/'+sample+'_barcode_noretrim_0_R2.fq'
    os.system('mv %s %s'%(fq1c, fq1))
    os.system('mv %s %s'%(fq2c, fq2))
    if config.RNA_PE==True:
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1},{2} -m {3}:{4} -a {5}  -o {6}/{7}_barcode_noretrim_1_R1.fq -p {6}/{7}_barcode_noretrim_1_R2.fq  {8} {9} 1>{6}/{7}_0.log""".format(nthreads,q5,q3,mlen1,mlen2,a1,trimmed_dir,sample,fq1,fq2), shell=True)
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1},{2} -m {3}:{4} -a {5}  -o {6}/{7}_barcode_R1.fq -p {6}/{7}_barcode_R2.fq  {6}/{7}_barcode_noretrim_1_R1.fq {6}/{7}_barcode_noretrim_1_R2.fq 1>{6}/{7}.log""".format(nthreads,q5,q3,mlen1,mlen2,a2,trimmed_dir,sample), shell=True)
    else:
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1} -m {2} -a {3}  -o {4}/{5}_barcode_noretrim_1_R2.fq {6} 1>{4}/{5}_0.log""".format(nthreads,q5,mlen2,a1,trimmed_dir,sample,fq2), shell=True)
        subprocess.call(CUTADAPT_PATH + """ -j {0} -q {1} -m {2} -a {3}  -o {4}/{5}_barcode_R2.fq {4}/{5}_barcode_noretrim_1_R2.fq 1>{4}/{5}.log""".format(nthreads,q5,mlen2,a2,trimmed_dir,sample), shell=True)
    return

def run_star(genome_dir, output_dir, sample, nthreads, PE):
    """ Align reads using STAR.
    """
    map_dir=output_dir+'/mapping'
    sample_dir=output_dir+'/split_sample'
    nthreads = int(nthreads)
    if PE==True:
        rc = subprocess.call(STAR_PATH + """ --genomeDir {0}/ --runThreadN {2} --readFilesIn {4}/{3}_barcode_R1.fq {4}/{3}_barcode_R2.fq --outFileNamePrefix {1}/{3}.""".format(genome_dir, map_dir, nthreads, sample, sample_dir), shell=True)
    else:
        rc = subprocess.call(STAR_PATH + """ --genomeDir {0}/ --runThreadN {2} --readFilesIn {4}/{3}_barcode_R2.fq --outFileNamePrefix {1}/{3}. --limitOutSJcollapsed 5000000""".format(genome_dir, map_dir, nthreads, sample, sample_dir), shell=True) # error encountered when --limitOutSJcollapsed is default 1000000
    # Add alignment stats to pipeline_stats
    with open('%s/%s.Log.final.out'%(map_dir,sample)) as f:
        for i in range(8):
            f.readline()
        unique_mapping = int(f.readline().split('\t')[1][:-1])
        for i in range(14):
            f.readline()
        multimapping = int(f.readline().split('\t')[1][:-1])
    with open( '%s/%s_pipeline_stats.txt'%(map_dir, sample), 'w') as f:
        #f.write('####### STAR mapping stats ########\n')
        f.write('uniquely_aligned\t%d\n' %unique_mapping)
        f.write('multimapping\t%d\n' %multimapping)
    return rc


def is_sorted_queryname(bamfile):
    """
    Check if bam fiel is sorted by read name.
    """
    samfile = pysam.AlignmentFile(bamfile, "rb")
    header= samfile.header
    samfile.close()
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                return True
    return False

def is_sorted_coordinate(bamfile):                                                                                                                             
    """
    Check if bam fiel is sorted by read name.
    """
    samfile = pysam.AlignmentFile(bamfile, "rb")
    header= samfile.header
    samfile.close()
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "coordinate"):
                return True
    return False


def sort_sam(output_dir, sample, nthreads, intype='sam'):
    """ Sort samfile by header (cell_barcodes, umi) 
    """
    map_dir=output_dir+'/mapping'
    nthreads = int(nthreads)
    rc = subprocess.call(SAMTOOLS_PATH + """ sort -n -@ {1} -T {0}/{2}.Aligned.sort -o {0}/{2}.Aligned.sorted.bam {0}/{2}.Aligned.out.{3}""".format(map_dir, nthreads, sample, intype), shell=True)
    os.remove("%s/%s.Aligned.out.%s"%(map_dir,sample, intype))
    return rc


def addTag(inbam, outbam):
    samfile = pysam.Samfile(inbam)
    new_samfile = pysam.Samfile(outbam,'wb', template=samfile)
    for read in samfile:
        cell_barcode = read.qname.split(':')[0]
        umi=read.qname.split(':')[1]
        read.tags+=[('CB',cell_barcode),('UB',umi)]
        new_samfile.write(read)
    samfile.close()
    new_samfile.close()