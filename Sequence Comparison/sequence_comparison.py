# COMP 561 - Final Project Code (Fall 2019)
# Dorian Desblancs - 260722712

# WARNING:  The program as presented takes about two hours to run
# ADVICE:   Go through main method and change the code to tailor it 
#           to your needs and/or understand what is going on

from Bio import SeqIO
import sys
import pandas as pd
import random
import numpy as np
import xlsxwriter

# String to list
def split(s):
    return [char for char in s]

# Extract all sequences from FASTA file
def fasta(f):
    record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    l = []
    for key, value in record_dict.items():
        sequence = record_dict.get(key)
        l_seq    = split(sequence)
        two_dim  = [key, l_seq]
        l.append(two_dim)
    return l

# Extract all sequences from XLSX file
def xlsx(f):
    sequences = []
    for index in range(1, 33, 1):
        sequence = pd.read_excel(f, skiprows=9, usecols=[index])
        seq_name = sequence._info_axis.values[0]
        sequence = sequence.values.tolist()
        flat = []
        for sublist in sequence:
            for item in sublist:
                item = item.replace(" ", "")
                flat.append(item)
        sequences.append([seq_name, flat])
    return sequences

# Following is adapted from Alexandr Levchuk's Needleman-Wunsch implementation
# His code can be found here:
# https://github.com/alevchuk/pairwise-alignment-in-python

def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

match_award      = 5
mismatch_penalty = -2
gap_penalty      = -30 # both for opening and extanding

def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty

"""
def finalize(seq1_name, seq2_name, align1, align2):
    align1 = align1[::-1]    #reverse sequence 1
    align2 = align2[::-1]    #reverse sequence 2
    
    i,j = 0,0
    
    #calcuate identity, score and aligned sequeces
    #symbol = ''
    #found = 0
    score = 0
    identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:                
            #symbol = symbol + align1[i]
            identity = identity + 1
            score += match_score(align1[i], align2[i])
    
        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += match_score(align1[i], align2[i])
            #symbol += ' '
            #found = 0
    
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':          
            #symbol += ' '
            score += gap_penalty
    
    identity = float(identity) / len(align1) * 100
    
    ov = round(identity, 3)
    sc = score

    return ov, sc
"""

def needle(seq1, seq2):
    seq1_name = seq1[0]
    seq2_name = seq2[0]
    seq1 = seq1[1]
    seq2 = seq2[1]

    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table
   
    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if (len(seq1[i-1]) == 1 and len(seq2[j-1]) == 1):
                match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1])
                delete = score[i - 1][j] + gap_penalty
                insert = score[i][j - 1] + gap_penalty
                score[i][j] = max(match, delete, insert)
            elif (len(seq1[i-1]) > 1 or len(seq2[j-1]) > 1):
                match = -np.inf
                delete = -np.inf
                insert = -np.inf
                score_param = -np.inf
                for l1 in range(len(seq1[i-1])):
                    for l2 in range(len(seq2[j-1])):
                        temp_match = score[i-1][j-1] + match_score(seq1[i-1][l1], seq2[j-1][l2])
                        temp_delete = score[i-1][j] + gap_penalty
                        temp_insert = score[i][j-1] + gap_penalty
                        temp_score = max(temp_match, temp_delete, temp_insert)
                        #if (temp_match > match):
                            #match = temp_match
                            #temp_nuc1 = seq1[i-1][l1]
                            #temp_nuc2 = seq2[j-1][l2]
                        if (temp_score > score_param):
                            score_param = temp_score
                score[i][j] = score_param
                #seq1[i-1] = temp_nuc1
                #seq2[j-1] = temp_nuc2
    """
    # Traceback and compute the alignment 
    align1, align2 = '', ''
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: # end toching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1]):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    ov, sc = finalize(seq1_name, seq2_name, align1, align2)
    
    return ov, sc
    """
    return score[m][n]

# main method: used to get command-line arguments and call all methods
def main(fa, xx, xl1, xl2):
    ll = fasta(fa)
    print("Finished gathering FASTA data...")
    seq= xlsx(xx)
    print("Finished gathering ancestral sequences...")
    print("Exporting Results to XLSX...")
    """
    xbook1 = xlsxwriter.Workbook(xl1)
    xsheet1 = xbook1.add_worksheet()
    xbook2 = xlsxwriter.Workbook(xl2)
    xsheet2 = xbook2.add_worksheet()
    column = 1
    for sequence in seq:
        row = 1
        print(sequence[0])
        xsheet1.write(0, column, int(sequence[0]))
        xsheet2.write(0, column, int(sequence[0]))
        for item in ll:
            ov, sc = needle(item, sequence)
            xsheet1.write(row, column, ov)
            xsheet2.write(row, column, sc)
            row = row + 1
        column = column + 1
    xbook1.close()
    xbook2.close()
    """
    xbook1 = xlsxwriter.Workbook(xl2)
    xsheet1 = xbook1.add_worksheet()
    column = 1
    for sequence in seq:
        row = 1
        print(sequence[0])
        xsheet1.write(0, column, int(sequence[0]))
        for item in ll:
            sc = needle(item, sequence)
            xsheet1.write(row, column, sc)
            row = row + 1
        column = column + 1
    xbook1.close()

main("Concatenated Fasta Files/COX1.fa", "Ancestral Sequences/cox1.xlsx", 'cox1_overlap.xlsx', 'cox1_score.xlsx')
print("Finished COX1 Results")
main("Concatenated Fasta Files/COX2.fa", "Ancestral Sequences/cox2.xlsx", 'cox2_overlap.xlsx', 'cox2_score.xlsx')
print("Finished COX2 Results")
main("Concatenated Fasta Files/COX3.fa", "Ancestral Sequences/cox3.xlsx", 'cox3_overlap.xlsx', 'cox3_score.xlsx')
print("Finished COX3 Results")