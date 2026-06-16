import math
from scipy.stats import binomtest
from operator import itemgetter
import time
import logomaker 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


residue_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
dipeptide_list = []
for residue1 in residue_list:
    for residue2 in residue_list:
        dipeptide = residue1 + residue2
        dipeptide_list.append(dipeptide)
tetrapeptide_list = []
for residue1 in residue_list:
    for residue2 in residue_list:
        for residue3 in residue_list:
            for residue4 in residue_list:
                tetrapeptide = residue1 + residue2 + residue3 + residue4
                tetrapeptide_list.append(tetrapeptide)

def fdr_correction(p_values, alpha=0.05):
    """
    Apply False Discovery Rate (FDR) correction to a list of p-values using the Benjamini-Hochberg procedure.
    
    Args:
        p_values (list): List of p-values.
        alpha (float): Desired significance level (default: 0.05).
        
    Returns:
        dictionary that applies correction to an input p-val
    """
    p_values_sorted = sorted(p_values)
    m = len(p_values_sorted)
    adjusted_p_values_sorted = [min(p * m / (i + 1), 1.0) for i, p in enumerate(p_values_sorted)]
    correction_dict = {}
    for pos in range(m):
        correction_dict[p_values_sorted[pos]] = adjusted_p_values_sorted[pos]
    adjusted_p_values = []
    for pos in range(m):
        adjusted_p_values.append(correction_dict[p_values[pos]])
    return correction_dict


        
def IC(residue_frequency_dict, pos):
    cumulative_value = 0
    for residue in residue_list:
        residue_frequency = residue_frequency_dict[residue][pos]
        if residue_frequency != 0.0:
            cumulative_value -= residue_frequency * math.log(residue_frequency, 2)
    IC_value = 4.32193 - cumulative_value
    return IC_value

def MI(residue_frequency_dict, dipeptide_frequency_dict, pos1, pos2):
    cumulative_value = 0.0
    for res1 in residue_list:
        for res2 in residue_list:
            res1_frequency = residue_frequency_dict[pos1][res1]
            res2_frequency = residue_frequency_dict[pos2][res2]
            dipeptide = res1 + res2
            dipeptide_frequency = dipeptide_frequency_dict[(pos1,pos2)][dipeptide] 
            if dipeptide_frequency != 0:
                cumulative_value += dipeptide_frequency * math.log(dipeptide_frequency / (res1_frequency * res2_frequency), 2)
    return cumulative_value

def format_tetrapeptide(tetrapeptide, pos1, pos2, pos3, pos4, peptide_length, res322):
    output_string = ''
    for pos in range(peptide_length):
        if pos == 8:
            output_string += res322
        elif pos == pos1 or pos == pos2 or pos == pos3 or pos == pos4:
            if pos == pos1:
                output_string += tetrapeptide[0]
            if pos == pos2:
                output_string += tetrapeptide[1]
            if pos == pos3:
                output_string += tetrapeptide[2]
            if pos == pos4:
                output_string += tetrapeptide[3]
            
        else:
            output_string += 'X'
    return output_string

def pattern_match(peptide, pattern):
    match = 1
    for pos in range(len(pattern)):
        if pattern[pos] != 'X' and peptide[pos] != 'X' and pattern[pos] != peptide[pos]:
            match = 0
    return match

#identifies all pairs of peptides sequences that differ by one residue
#comparing each peptide against all other peptides gets very slow as number of peptides increases
#this function speeds things up dramatically by only comparing peptides where either the first half or second half match
#all pairs of peptides with exacly one mismatch will have match at either the N-terminal or C-terminal ends 
def find_neighbor_peptides_by_dict(peptide_list):
    peptide_neighbor_dict = {}
    for peptide in peptide_list:
        peptide_neighbor_dict[peptide] = []
    midpoint = int(len(peptide_list[0]) / 2)
    N_pep_dict = {} #dictionary of sequences of the N-terminal half of the peptide
    C_pep_dict = {} #same as above, except for C-terminal half
    for pep_pos in range(len(peptide_list)):
        peptide = peptide_list[pep_pos]
        N_pep = peptide[:midpoint]
        C_pep = peptide[midpoint:]
        if N_pep not in N_pep_dict:
            N_pep_dict[N_pep] = [peptide]
        else:
            N_pep_dict[N_pep].append(peptide)
        if C_pep not in C_pep_dict:
            C_pep_dict[C_pep] = [peptide]
        else:
            C_pep_dict[C_pep].append(peptide)
    comparison_dict = {} #keeps track of comparisons made for N-terminal matches to avoid double-counting anything; might not be needed
    for N_pep in list(N_pep_dict.keys()):
        for pep1_pos in range(len(N_pep_dict[N_pep])):
            for pep2_pos in range(len(N_pep_dict[N_pep])):
                if pep2_pos > pep1_pos:
                    peptide1 = N_pep_dict[N_pep][pep1_pos]
                    peptide2 = N_pep_dict[N_pep][pep2_pos]
                    mismatches = 0
                    for pos in range(len(peptide1)):
                        if peptide1[pos] != peptide2[pos]:
                            mismatches += 1
                        if mismatches > 1:
                            break
                    if mismatches == 1:
                        #print(peptide1, peptide2)
                        peptide_neighbor_dict[peptide1].append(peptide2)
                        peptide_neighbor_dict[peptide2].append(peptide1)
                        comparison_dict[(peptide1, peptide2)] = 1
                        comparison_dict[(peptide2, peptide1)] = 1
    for C_pep in list(C_pep_dict.keys()):
        for pep1_pos in range(len(C_pep_dict[C_pep])):
            for pep2_pos in range(len(C_pep_dict[C_pep])):
                if pep2_pos > pep1_pos:
                    peptide1 = C_pep_dict[C_pep][pep1_pos]
                    peptide2 = C_pep_dict[C_pep][pep2_pos]
                    if (peptide1, peptide2) not in comparison_dict:
                        mismatches = 0
                        for pos in range(len(peptide1)):
                            if peptide1[pos] != peptide2[pos]:
                                mismatches += 1
                            if mismatches > 1:
                                break
                        if mismatches == 1:
                            #print(peptide1, peptide2)
                            peptide_neighbor_dict[peptide1].append(peptide2)
                            peptide_neighbor_dict[peptide2].append(peptide1)
    return peptide_neighbor_dict
                  
def pattern_create(seq1, seq2):
    pattern = ''
    for pos in range(len(seq1)):
        if seq1[pos] == seq2[pos]:
            pattern += seq1[pos]
        else:
            pattern += 'X'
    return pattern

def pattern_compare(pattern1, pattern2):
    matches = 0
    mismatches = 0
    for pos in range(len(pattern1)):
        if pattern1[pos] != 'X' and pattern2[pos] != 'X':
            if pattern1[pos] == pattern2[pos]:
                matches += 1
            else:
                mismatches += 1
    return (matches, mismatches)


if len(sys.argv) < 2:
    print('please give the input filename as a command line argument')
else:
    data_filename = sys.argv[1]
    truncated_filename = data_filename[:-4]
    pattern_filename = truncated_filename + '_4res_patterns.txt' 
    pattern_file = open(pattern_filename, 'w')
    plot_peptide_list = []
    
    for res322 in ['G', 'A', 'R', 'P']:
        #print('residue at position 322 = %s' %res322)
        read_threshold = 1
        tetrapeptide_count_threshold = 3
        total_peptides = 0
        peptide_list = []
        data_file = open(data_filename, 'r')
        for raw_line in data_file:
            line = raw_line.strip()
            if line:
                DNA, raw_DNA_count, peptide = line.split('\t')
                DNA_count = int(raw_DNA_count)
                if DNA_count >= read_threshold:
                    total_peptides += 1
                    if peptide[8] == res322:
                        peptide_list.append(peptide)

            
        peptide_dict = {}
        dipeptide_dict = {}
        tetrapeptide_dict = {}
        peptide_length = len(peptide_list[0])
        if peptide_length == 6:
            exclude_pos_list = []
        if peptide_length == 7:
            exclude_pos_list = [4]
            peptide_offset = 231
        if peptide_length == 9:
            exclude_pos_list = [1, 2, 6]
        if peptide_length == 12:
            exclude_pos_list = [1, 3, 5, 6, 8, 10] #don't look for patterns involving residue322 which is only allowed to be G, A, R, or P
            peptide_offset = 314


        for pos in range(peptide_length):
            peptide_dict[pos] = {}
            for residue in residue_list:
                peptide_dict[pos][residue] = 0
        for pos1 in range(peptide_length):
            for pos2 in range(peptide_length):
                if pos2 > pos1:
                    dipeptide_dict[(pos1,pos2)] = {}
                    for dipeptide in dipeptide_list:
                        dipeptide_dict[(pos1,pos2)][dipeptide] = 0
        for pos1 in range(peptide_length):
            for pos2 in range(peptide_length):
                for pos3 in range(peptide_length):
                    for pos4 in range(peptide_length):
                        if pos2 > pos1 and pos3 > pos2 and pos4 > pos3 and pos1 not in exclude_pos_list and pos2 not in exclude_pos_list and pos3 not in exclude_pos_list and pos4 not in exclude_pos_list:
                            tetrapeptide_dict[(pos1,pos2,pos3,pos4)] = {}
                            for tetrapeptide in tetrapeptide_list:
                                tetrapeptide_dict[(pos1,pos2,pos3,pos4)][tetrapeptide] = 0




        for peptide in peptide_list:
            for pos in range(peptide_length):
                peptide_dict[pos][peptide[pos]] += 1
            for pos1 in range(peptide_length):
                for pos2 in range(peptide_length):
                    if pos2 > pos1:
                        dipeptide_dict[(pos1, pos2)][peptide[pos1] + peptide[pos2]] += 1
        for peptide in peptide_list:
            for pos1 in range(peptide_length):
                for pos2 in range(peptide_length):
                    for pos3 in range(peptide_length):
                        for pos4 in range(peptide_length):
                            if pos2 > pos1 and pos3 > pos2 and pos4 > pos3 and pos1 not in exclude_pos_list and pos2 not in exclude_pos_list and pos3 not in exclude_pos_list and pos4 not in exclude_pos_list:
                                tetrapeptide = peptide[pos1] + peptide[pos2] + peptide[pos3] + peptide[pos4]
                                tetrapeptide_dict[(pos1,pos2,pos3,pos4)][tetrapeptide] += 1
        ##for residue in residue_list:
        ##    output_line = '%s\t' %residue
        ##    for pos in range(peptide_length):
        ##        output_line += '%5.2f\t' %(float(peptide_dict[pos][residue]) / float(total_peptides))
        ##    print(output_line)
        ##print()

        print('%d of %d peptides have %s at position 322' %(len(peptide_list), total_peptides, res322))

        residue_frequency_dict = {}
        for pos in range(peptide_length):
            residue_frequency_dict[pos] = {}
            for res in residue_list:
                residue_frequency_dict[pos][res] = float(peptide_dict[pos][res]) / float(total_peptides)
        ##for pos in range(peptide_length):
        ##    IC_for_pos = IC(residue_frequency_dict, pos)
        ##    print('%d\t%5.2f' %(pos, IC_for_pos))
        ##print(IC(residue_frequency_dict, pos))

        dipeptide_frequency_dict = {}
        for pos1 in range(peptide_length):
            for pos2 in range(peptide_length):
                if pos2 > pos1:
                    dipeptide_frequency_dict[(pos1,pos2)] = {}
                    for dipeptide in dipeptide_list:
                        dipeptide_frequency_dict[(pos1,pos2)][dipeptide] = float(dipeptide_dict[(pos1,pos2)][dipeptide]) /  float(total_peptides)

        ##for pos1 in range(peptide_length):
        ##    for pos2 in range(peptide_length):
        ##        if pos2 > pos1:
        ##            for res1 in residue_list:
        ##                for res2 in residue_list:
        ##                    dipeptide = res1 + res2
        ##                    if dipeptide_dict[(pos1, pos2)][dipeptide] >= 3.0:
        ##                        enrichment = dipeptide_frequency_dict[(pos1,pos2)][dipeptide] / (residue_frequency_dict[pos1][res1] * residue_frequency_dict[pos2][res2])
        ##                        if enrichment >= 2.0:
        ##                            print('%d\t%s\t%d\t%s\t%5.2f' %(pos1, res1, pos2, res2, enrichment))
        enriched_tetramer_dict = {}
        comparisons = 0
        for pos1 in range(peptide_length):
            for pos2 in range(peptide_length):
                for pos3 in range(peptide_length):
                    for pos4 in range(peptide_length):
                        if pos2 > pos1 and pos3 > pos2 and pos4 > pos3 and pos1 not in exclude_pos_list and pos2 not in exclude_pos_list and pos3 not in exclude_pos_list and pos4 not in exclude_pos_list:
                            for res1 in residue_list:
                                for res2 in residue_list:
                                    for res3 in residue_list:
                                        for res4 in residue_list:
                                            tetrapeptide = res1 + res2 + res3 + res4
                                            if tetrapeptide_dict[(pos1, pos2, pos3, pos4)][tetrapeptide] >= 1:
                                                formatted_tetrapeptide = format_tetrapeptide(tetrapeptide, pos1, pos2, pos3, pos4, peptide_length, res322)
                                                observed_tetrapeptide = tetrapeptide_dict[(pos1, pos2, pos3, pos4)][tetrapeptide]
                                                expected_tetrapeptide = float(total_peptides) * (residue_frequency_dict[pos1][res1] * residue_frequency_dict[pos2][res2] * residue_frequency_dict[pos3][res3] * residue_frequency_dict[pos4][res4])
                                                expected_fraction = (residue_frequency_dict[pos1][res1] * residue_frequency_dict[pos2][res2] * residue_frequency_dict[pos3][res3] * residue_frequency_dict[pos4][res4])
                                                enrichment = float(observed_tetrapeptide) / expected_tetrapeptide
                                                if enrichment > 1:
                                                    dipeptide12 = res1 + res2
                                                    dipeptide34 = res3 + res4
                                                    a = tetrapeptide_dict[(pos1, pos2, pos3, pos4)][tetrapeptide]
                                                    b = dipeptide_dict[(pos1,pos2)][dipeptide12] - a
                                                    c = dipeptide_dict[(pos3,pos4)][dipeptide34] - a
                                                    d = total_peptides - a - b -c
                                                    table = np.array([[a, b], [c, d]])
                                                    #result = fisher_exact(table, alternative='greater')
                                                    result = binomtest(observed_tetrapeptide, total_peptides, expected_fraction, alternative='greater') 
                                                    pvalue = result.pvalue
                                                    #print(pvalue)
                                                    enriched_tetramer_dict[formatted_tetrapeptide] = (pvalue, enrichment, tetrapeptide_dict[(pos1, pos2, pos3, pos4)][tetrapeptide])
        sorted_enrichment_list = sorted(enriched_tetramer_dict.items(), key=lambda x:x[1], reverse=False)
        pval_list = []
        for item in sorted_enrichment_list:
            pval_list.append(item[1][0])
            #print(item[1][0])
        pval_correction_dict = fdr_correction(pval_list)
        pattern_file.write('#%d of %d peptides have %s at position 322\n' %(len(peptide_list), total_peptides, res322))
        filtered_pattern_list = sorted_enrichment_list[:100]
        final_pattern_list = []
        for item in filtered_pattern_list:
            pval = item[1][0]
            #print(pval)
            corrected_pval = pval_correction_dict[pval]
            #print(corrected_pval)
            if corrected_pval < 0.05 and item[1][2] >= tetrapeptide_count_threshold:
                print('%s\t%6.2e\t%5.2f\t%d' %(item[0], corrected_pval, item[1][1], item[1][2]))
                pattern_file.write('%s\t%s\t%6.2e\t%5.2f\t%d\n' %(truncated_filename, item[0], corrected_pval, item[1][1], item[1][2]))
                final_pattern_list.append(item)

        data_file.close()

        #this section gives a score to each peptide for how well it matches the enriched patterns.
        #we use this score to help decide which peptide sequences to build and experimentally validate
        peptide_pattern_score_dict = {}
        for peptide in peptide_list:
            pattern_match_score = 0
            for pattern in final_pattern_list:
                if pattern_match(peptide, pattern[0]) == 1:
                    if pattern[1][0] > 0:
                        pattern_match_score += -1.0 * math.log(pattern[1][0], 10)
                    else:
                        pattern_match_score += 312.0 # p values below abount 10^-312 are set to 0 and cause math.log to throw and error
            peptide_pattern_score_dict[peptide] = pattern_match_score

        sorted_peptide_score_list = sorted(peptide_pattern_score_dict.items(), key=lambda x:x[1], reverse=True)
        residue_pos_dict = {}
        residue_pos_prob_dict = {}
        IC_scaled_residue_pos_prob_dict = {}
        for res in residue_list:
            residue_pos_dict[res] = {}
            residue_pos_prob_dict[res] = {}
            IC_scaled_residue_pos_prob_dict[res] = {}
            for pos in range(peptide_length):
                residue_pos_dict[res][pos] = 0
                residue_pos_prob_dict[res][pos + peptide_offset] = 0.0
                IC_scaled_residue_pos_prob_dict[res][pos + peptide_offset] = 0.0

        output_filename = truncated_filename + '_top_scores_322%s.txt' %res322
        output_file = open(output_filename, 'w')
        
        for item in sorted_peptide_score_list:#[:100]:
            #print('%s\t%5.2f' %(item[0], item[1]))
            if item[1] > 0:
                output_file.write('%s\t%5.2f\n' %(item[0], item[1]))
        output_file.close()
