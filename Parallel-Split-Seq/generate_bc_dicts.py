# make barcode dictionaries
import os
import subprocess
import sys

import pandas as pd
from collections import defaultdict
import gzip
from numpy import unique
import numpy as np
import pickle

bases=['A','G','C','T']

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

def free_divergence_matrix(s1, s2):
    assert len(s1) == len(s2), 'Free divergence requires strings of equal length.'
    n = len(s1) + 1
    M = np.zeros((n, n))

    # First row and column
    M[0, :] = np.arange(n)
    M[:, 0] = np.arange(n)

    # Inner square
    for i in range(1, n-1):
        for j in range(1, n-1):
            if s1[i-1] == s2[j-1]:
                diag_penalty = M[i-1, j-1]
            else:
                diag_penalty = M[i-1, j-1] + 1
            M[i, j] = min(diag_penalty,
                          M[i-1, j] + 1,
                          M[i, j-1] + 1)

    # Last row
    for i in range(1, n-1):
        if s1[n-2] == s2[i-1]:
            diag_penalty = M[n-2, i-1]
        else:
            diag_penalty = M[n-2, i-1] + 1
        M[n-1, i] = min(diag_penalty,
                        M[n-2, i] + 1,
                        M[n-1, i-1])  # No penalty along last row

    # Last column
    for i in range(1, n-1):
        if s1[i-1] == s2[n-2]:
            diag_penalty = M[i-1, n-2]
        else:
            diag_penalty = M[i-1, n-2] + 1
        M[i, n-1] = min(diag_penalty,
                        M[i-1, n-1],  # No penalty along last col
                        M[i, n-2] + 1)

    # Last elt
    if s1[n-2] == s2[n-2]:
        diag_penalty = M[n-2, n-2]
    else:
        diag_penalty = M[n-2, n-2] + 1
    M[n-1, n-1] = min(diag_penalty,
                      M[n-2, n-1],  # No penalty along last col
                      M[n-1, n-2])  # No penalty along last row
    return M

def free_divergence(s1, s2):
    M = free_divergence_matrix(s1, s2)
    return M[-1, -1]

def generate_edit_dict(barcode_dict, mer, method):
    edit_dict = {0:defaultdict(list),1:defaultdict(list),2:defaultdict(list),3:defaultdict(list)}
    for bc in barcode_dict:
        for cur_mer in mer:
            if method=='Levenshtein':
                lev_dist = levenshteinDistance(bc,cur_mer)
            elif method=='FreeDivergence':
                lev_dist = free_divergence(bc,cur_mer)
            if lev_dist==0:
                edit_dict[0][cur_mer].append(bc)
            elif lev_dist==1:
                edit_dict[1][cur_mer].append(bc)
            elif lev_dist==2:
                edit_dict[2][cur_mer].append(bc)
            elif lev_dist==3:
                edit_dict[3][cur_mer].append(bc)
        # print(bc,end=' ')
    return edit_dict

def generate_edit_dict_mix(barcode_dict, mer, method):
    edit_dict = {0:defaultdict(list),1:defaultdict(list),2:defaultdict(list)}
    for bc in barcode_dict:
        bcl = len(bc)
        for cur_mer in mer[bcl]:
            if method=='Levenshtein':
                lev_dist = levenshteinDistance(bc,cur_mer)
            elif method=='FreeDivergence':
                lev_dist = free_divergence(bc,cur_mer)
            if lev_dist==0:
                edit_dict[0][cur_mer].append(bc)
            elif lev_dist==1:
                edit_dict[1][cur_mer].append(bc)
            elif lev_dist==2:
                edit_dict[2][cur_mer].append(bc)
            #elif lev_dist==3:
            #    edit_dict[3][cur_mer].append(bc)
        # print(bc,end=' ')
    return edit_dict

def write_pkl(method, bcset):
    '''
    #mer6 = convert_degen_seq_to_list('NNNNNN')
    mer8 = convert_degen_seq_to_list('NNNNNNNN')
    mer9 = convert_degen_seq_to_list('NNNNNNNNN')
    mer10 = convert_degen_seq_to_list('NNNNNNNNNN')
    mer11 = convert_degen_seq_to_list('NNNNNNNNNNN')
    
    bc1_8nt = pd.read_csv(bcset+'/barcode1_part.txt',names=['barcode'],index_col=0,sep='\t').barcode
    edit_dict_1 = generate_edit_dict(bc1_8nt, mer8, method)
    with open(bcset+'/bc_dict_1_part_%s.pkl'%method, 'wb') as f:
        pickle.dump(edit_dict_1, f)
    
    mer = {8:mer8,9:mer9,10:mer10,11:mer11}
    bc1_mix = pd.read_csv(bcset+'/barcode1.txt',names=['barcode'],index_col=0,sep='\t').barcode
    edit_dict_1 = generate_edit_dict_mix(bc1_mix, mer, method)
    with open(bcset+'/bc_dict_1_%s.pkl'%method, 'wb') as f:
        pickle.dump(edit_dict_1, f)
    '''
    
    bc1_8nt = pd.read_csv(bcset+'/barcode1.txt',names=['barcode'],index_col=0,sep='\t').barcode
    edit_dict_1 = generate_edit_dict(bc1_8nt, mer8, method)
    with open(bcset+'/bc_dict_1_%s.pkl'%method, 'wb') as f:
        pickle.dump(edit_dict_1, f)
    bc2_8nt = pd.read_csv(bcset+'/barcode2.txt',names=['barcode'],index_col=0,sep='\t').barcode
    edit_dict_2 = generate_edit_dict(bc2_8nt, mer8, method)
    with open(bcset+'/bc_dict_2_%s.pkl'%method, 'wb') as f:
        pickle.dump(edit_dict_2, f)
    
    bc3_6nt = pd.read_csv(bcset+'/barcode3.txt',names=['barcode'],index_col=0,sep='\t').barcode
    edit_dict_3 = generate_edit_dict(bc3_6nt, mer6, method)
    with open(bcset+'/bc_dict_3_%s.pkl'%method, 'wb') as f:
        pickle.dump(edit_dict_3, f)
        
    
if __name__ == '__main__':
    import sys
    bcset=sys.argv[1]
    write_pkl('FreeDivergence',bcset)
