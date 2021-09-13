# Converting bam to read information
# Cell barcode, UMI, gene alignment, and counts
import os
import sys
from multiprocessing import Process

import pandas as pd
from collections import defaultdict
import gzip
#from pylab import *
import numpy as np
import pickle

#import HTSeq
import pysam

PATH = os.path.dirname(__file__)

bases = list('ACGT')
def single_mut_seqs(seq,n):
    """ Return a list of all sequences with a 1 nt mutation (no ins/del)."""
    mut_seqs = []
    for i in range(n):
        for b in bases:
            if b!=seq[i]:
                mut_seqs.append(seq[:i] + b + seq[i+1:])
    return mut_seqs

def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""

    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
        if diffs>1:
            break
    return diffs

def collapse_umis(kmer_list,counts):
    """ Collapses UMIs that are within 1 hamming dist of each other. Collapses by adding together counts
    of both similar UMIs."""
    kmer_len = len(kmer_list[0])
    kmer_seqs_dict = dict(zip(kmer_list,np.zeros(len(kmer_list))))
    ham_dict = defaultdict(list)
    seq_len = len(kmer_list)
    if seq_len > 50:
        for seq in kmer_list:
            mut_seqs = single_mut_seqs(seq,kmer_len)
            for ms in mut_seqs:
                try:
                    kmer_seqs_dict[ms]
                except:
                    pass
                else:
                    ham_dict[ms].append(seq)
    else:
        for i in range(seq_len):
            for j in range(i+1,seq_len):
                s1,s2 = kmer_list[i],kmer_list[j]
                ham_dist = hamdist(s1,s2)
                if ham_dist <= 1:
                    ham_dict[s1].append(s2)
                    ham_dict[s2].append(s1)
    kmer_counts = dict(zip(kmer_list,counts))
    valid_kmers = []
    for kmer in kmer_list:
        cur_count = kmer_counts[kmer]
        ham_matches = ham_dict[kmer]
        found_match = False
        for hm in ham_matches:
            try:
                match_count = kmer_counts[hm]
                if match_count > (cur_count):
                    found_match = True
                    kmer_counts[hm] += cur_count
                    del kmer_counts[kmer]
                    break
            except KeyError:
                pass
        #if not found_match:
        #    valid_kmers.append(kmer)    
    return kmer_counts

def collapse_umis_dataframe(df):
    """ Collapses similar (<1 hamming dist) UMIs with the same cell_barcode-gene combination """
    counts = df.groupby(['gene','umi']).size().reset_index()
    counts.columns = ['gene','umi','counts']
    counts_df = counts.groupby('gene').apply(lambda x: pd.Series(collapse_umis(list(x['umi']),list(x['counts']))))
    if not(type(counts_df) == pd.core.series.Series):    # Need to do this in case counts_df only has one
        umi = counts_df.columns[0]
        umi_counts = counts_df.iloc[:,0]
        counts_df['umi']=umi
        counts_df['counts'] = umi_counts
        counts_df = counts_df[['umi','counts']]
    return counts_df

def split_bam(output_dir, sample, nthreads):
    samfile = pysam.Samfile(output_dir + '/mapping/'+sample+'.Aligned.sorted.bam')
    samfile_chunks = {}
    for i in range(nthreads):
        samfile_chunks[i] = pysam.Samfile(output_dir + '/mapping/'+sample+'.Aligned.sorted.chunk%d.bam' %(i+1),'wb', template=samfile)
        
    # Get the total number of aligned reads:
    aligned_reads = 0
    with open(output_dir +'/mapping/'+sample+'_pipeline_stats.txt') as f:   
        aligned_reads += int(f.readline()[:-1].split('\t')[1])
        aligned_reads += int(f.readline()[:-1].split('\t')[1])

    # Number of reads per file. Round up to ensure we don't write (nthreads +1) files.
    reads_per_chunk = int(np.ceil(aligned_reads/nthreads))*2
    print(aligned_reads, reads_per_chunk)
    c = 0
    prev_cell_barcode = ''
    d = 0
    reads = []
    for read in samfile:
        if not read.is_secondary:
            reads.append(read)
            cell_barcode = read.qname.split(':')[0]
            # Only write reads once all reads for a cell have been loaded to avoid
            # splitting one transcriptome into multiple files:
            if cell_barcode!=prev_cell_barcode:
                chunk = int(np.floor(c/reads_per_chunk))
                for r in reads[:-1]:
                    samfile_chunks[chunk].write(r)
                reads = reads[-1:]
                d += 1
            prev_cell_barcode = cell_barcode
            c+=1
            
    for i in range(nthreads):
        samfile_chunks[i].close()
    samfile.close()

def molecule_info_chunk(transcriptome_dir, output_dir, sample, config, chunk=None, gtf_dict_stepsize=10000):
    """ Gets the molecular info for each UMI in a bamfile.
    """


    bamfile = output_dir + '/mapping/'+sample+ '.Aligned.sorted.bam'

    if chunk is None:
        bamfile = output_dir + '/mapping/'+sample+ '.Aligned.sorted.bam'
        output_filename = output_dir +'/molecule_info/'+sample+ '_read_assignment.csv'
        outstats = output_dir +'/molecule_info/'+sample+ '_bam_stats.txt'
    else:
        bamfile = output_dir +'/mapping/'+sample+ '.Aligned.sorted.chunk%d.bam' %chunk
        output_filename = output_dir +'/molecule_info/'+sample+ '_read_assignment.chunk%d.csv' %chunk
        outstats = output_dir +'/molecule_info/'+sample+ '_bam_stats.chunk%d.txt'%chunk

    with open(transcriptome_dir +'/gene_info.pkl', 'rb') as f:
        gene_info = pickle.load(f)

    gene_dict = gene_info['gene_bins']
    exon_gene_start_end_dict = gene_info['genes_to_exons']
    start_dict = gene_info['gene_starts']
    end_dict = gene_info['gene_ends']
    gene_id_to_name = gene_info['gene_id_to_name']
    gene_id_to_genome = gene_info['gene_id_to_genome']

    gene_id_to_name['Intergenic_mm']='Intergenic'
    gene_id_to_name['Ambiguous_mm']='Ambiguous'
    gene_id_to_genome['Intergenic_mm']='mm'
    gene_id_to_genome['Ambiguous_mm']='mm'
    
    gene_id_to_name['Intergenic_hg']='Intergenic'
    gene_id_to_name['Ambiguous_hg']='Ambiguous'
    gene_id_to_genome['Intergenic_hg']='hg'
    gene_id_to_genome['Ambiguous_hg']='hg'
    
    gene_id_to_name['Intergenic_dm']='Intergenic'
    gene_id_to_name['Ambiguous_dm']='Ambiguous'
    gene_id_to_genome['Intergenic_dm']='dm'
    gene_id_to_genome['Ambiguous_dm']='dm'
    
    

    samfile = pysam.Samfile(bamfile)
    chrom_dict = pd.DataFrame(samfile.header['SQ']).SN.to_dict()
    samfile.close()

    def get_gene(read):
        """ Returns a list of genes that read alignment maps to.
        The current implementation is stranded and will not return a gene on the opposite
        strand.
        Input: a read object from pysam Alignment file
        Output: a list of genes (or an empty list if read doesn't map to any gene)
        """
        if read.is_reverse:
            strand = '-'
        else:
            strand = '+'

        # Assign each read into a "bin" by rounding the position (<chrom>:<pos>:<strand>) to the nearest
        # gtf_dict_stepsize (default of 10kb).
        read_pos = chrom_dict[read.tid]+':'+str(int(np.floor(read.positions[0]/gtf_dict_stepsize)*gtf_dict_stepsize))+':'+strand

        # The gene_dict contains all genes that overlap a given "bin". This does not guarantee that the read actually overlaps
        # these genes, but reduces the list of possible genes to a few from >20k.
        potential_genes = gene_dict[read_pos]
        
        
        # Now we explicity check to see if the read maps to gene(s) from potential_genes.
        # The current implementation is strict such the full read must be contained in the gene and 
        # no part of the read can occur outside the gene.
        matching_genes = []
        for g in potential_genes:
            gene_start,gene_end = start_dict[g],end_dict[g]
            if (gene_start<=read.reference_start) and ((read.reference_start+abs(read.tlen))<=gene_end):
                matching_genes.append(g)
        return matching_genes

    def get_genome(read):
        chr_name=read.reference_id
        genome=chrom_dict[chr_name].split('_')[0]
        return genome

    def check_exon_alignment(read,gene):
        """ Returns the fraction of bases in a read that align to at least one exonic base
        """
        # TODO: Write a faster implementation
        possible_exons = exon_gene_start_end_dict[gene]
        k = len(possible_exons)

        # pysam is zero based and gtf is one based
        aligned_positions = np.arange(read.reference_start,read.reference_start+abs(read.tlen))

        align_array = np.zeros(len(aligned_positions),dtype=bool)

        for start,end in possible_exons.items():
            #print(read.positions)
            #print(start,end)
            align_array = align_array | ((aligned_positions>=(start)) & (aligned_positions<=(end)))
        return np.mean(align_array)

    # Collapse similar UMIs for the same gene-cell_barcode combination. Write output info for
    # each UMI (cell_barcode, gene, UMI, count) to a file (read_assignment.csv).
    with open(output_filename,'w') as f:
        f.write('cell_barcode,genome,gene,gene_name,umi,counts,exonic\n') 


    min_mapq=config.min_mapq
    max_flen=config.max_flen


    samfile = pysam.Samfile(bamfile)
    c = 0
    cell_barcodes = []
    #genomes = []
    genes = []
    umis = []
    exonic_assignments = []
    NUM_READS = 0
    NUM_UNIQ = 0
    NUM_CHRM=0
    NUM_PAIRED=0
    NUM_USABLE=0
    next_cell = False
    species_counts = {}
    for read in samfile:
        rname  = str(read.reference_name)
        NUM_READS += 1
        if (not read.is_secondary):
            NUM_UNIQ+=1
            if ("chrM" in rname) or ("mito" in rname):
                NUM_CHRM+=1
                continue
            if read.is_unmapped or read.mate_is_unmapped: continue
            if read.is_proper_pair and (not read.is_reverse):
                NUM_PAIRED+=1
                if (abs(read.isize) <= max_flen) and (read.mapq >= min_mapq) :
                    NUM_USABLE += 1
                    genome = get_genome(read)
                    #genomes.append(genome)
                    gene_match = get_gene(read)    
                    if not gene_match:
                        exonic_assignments.append(False)
                        genes.append('Intergenic_'+genome)
                    elif len(gene_match)==1:
                        exonic_assignments.append(check_exon_alignment(read,gene_match[0])>0.5)
                        genes.append(gene_match[0])
                    else:
                        exonic_assignments.append(False)
                        genes.append('Ambiguous_'+genome)
                    seqname=read.qname
                    cell_barcodes.append(seqname.split(':')[0].split('_')[-1])
                    umis.append(seqname.split(':')[1])

                    c+=1

                    # Only want to collapse UMIs and write reads to file once, the 
                    # all the reads from a cell have been loaded. As soon as we see 
                    # a read from a different cell than the previous read, we process
                    # all the reads from the first cell.
                    if c>1:
                            if (cell_barcodes[-1]!=cell_barcodes[-2]):
                                next_cell=True

                        # Process reads from one cell: collapse UMIs and write to file
                    if next_cell:
                            df = pd.DataFrame({'cell_barcodes':cell_barcodes[:-1],
                                               'umi':umis[:-1],
                                               'gene':pd.Series(genes[:-1]),
                                               'exonic':exonic_assignments[:-1]})
                            df_collapsed = collapse_umis_dataframe(df)
                            df_collapsed = df_collapsed.reset_index()
                            df_collapsed.columns = ['gene','umi','counts']
                            exonic_fraction = df.groupby(['gene','umi']).exonic.mean()
                            df_collapsed.loc[:,'exonic'] = exonic_fraction.loc[list(zip(df_collapsed.gene.values,
                                                                                     df_collapsed.umi.values))].values>0.5
                            df_collapsed['cell_barcode'] = cell_barcodes[0]
                            df_collapsed.loc[:,'gene_name'] = df_collapsed.gene.apply(lambda s:gene_id_to_name[s])
                            df_collapsed.loc[:,'genome'] = df_collapsed.gene.apply(lambda s:gene_id_to_genome[s])

                            df_collapsed[['cell_barcode','genome','gene','gene_name','umi','counts','exonic']].to_csv(output_filename,
                                                                           header=False,
                                                                           index=False,
                                                                           mode='a')

                            # Reset lists, but keep the first read from the new cell
                            cell_barcodes = cell_barcodes[-1:]
                            genes = genes[-1:]
                            umis = umis[-1:]
                            exonic_assignments = exonic_assignments[-1:]
                            c = 1
                            next_cell=False

        
        if NUM_READS%100000==0:
            print('Processed %d aligned reads...' %NUM_READS)
            sys.stdout.flush()
    samfile.close()
    with open(outstats,'w') as output:
        output.write("number of totol reads\t%d\n"%NUM_READS)
        output.write("number of uniquely mapped reads\t%d\n"%NUM_UNIQ)
        output.write("number of chrM reads\t%d\n"%NUM_CHRM)
        output.write("number of proper pairs\t%d\n"%NUM_PAIRED)
        output.write("number of usable pairs\t%d"%NUM_USABLE)


def join_read_assignment_files(output_dir, sample, nthreads):
    filenames = [output_dir +'/molecule_info/'+sample+ '_read_assignment.chunk%d.csv' %i for i in range(1,nthreads+1)]
    with open(output_dir + '/molecule_info/'+sample+ '_read_assignment.csv', 'w') as outfile:
        # Write header
        outfile.write('cell_barcode,genome,gene,gene_name,umi,counts,exonic\n')
        for fname in filenames:
            with open(fname) as infile:
                # Don't copy header line from each file:
                infile.readline()
                for line in infile:
                    outfile.write(line)

def join_atac_stats(output_dir,sample,nthreads):
    filenames = [output_dir +'/molecule_info/'+sample+ '_bam_stats.chunk%d.txt' %i for i in range(1,nthreads+1)]
    stats={}
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                key,value=line.strip().split('\t')
                stats.setdefault(key, 0)
                stats[key]+=int(value) 
    with open(output_dir +'/molecule_info/'+sample+ '_bam_stats.txt','w') as outfile:
        for key in stats.keys():
            outfile.write(key+'\t%d\n'%stats[key])

def molecule_info(transcritome_dir, output_dir, sample, config, nthreads):
    """ Gets molecular info for a bam file. Splits the bamfile into 
    nthread chunks and runs in parallel """
    nthreads = int(nthreads)

    split_bam(output_dir, sample, nthreads)

    Pros = []
    for i in range(1,nthreads+1):
        print('Starting thread %d' %i)
        p = Process(target=molecule_info_chunk, args=(transcritome_dir, output_dir,sample, config, i))
        Pros.append(p)
        p.start()
    for t in Pros:
        t.join()

    join_read_assignment_files(output_dir, sample, nthreads)
    join_atac_stats(output_dir,sample,nthreads)

    # Log the total reads mapped to transcriptome and total UMIs
    total_read_count = 0
    mapped_to_gene = 0
    with open(output_dir +'/molecule_info/'+sample+ '_read_assignment.csv') as f:
        for i, l in enumerate(f):
            if i>0:
                total_read_count += int(l.split(',')[-2])
                if l.split(',')[3]!='Intergenic' and l.split(',')[3]!='Ambiguous':
                    mapped_to_gene+=int(l.split(',')[-2])

    with open(output_dir +'/molecule_info/'+sample+ '_pipeline_stats.txt', 'w') as f:
        f.write('total_read_count\t%d\n' %total_read_count)
        f.write('mapped_to_gene\t%d\n' %mapped_to_gene)
        f.write('total_umis\t%d' %i)

    # Remove temporary split files:
    
    for i in range(1,int(nthreads)+1):
        os.remove(output_dir +'/molecule_info/'+sample+ '_read_assignment.chunk%d.csv' %i)
        os.remove(output_dir +'/mapping/'+sample+ '.Aligned.sorted.chunk%d.bam' %i)
        os.remove(output_dir +'/molecule_info/'+sample+ '_bam_stats.chunk%d.txt' %i)
    
