#this script identifies enriched 4 residue motifs in selected Bxb1 helix sequences
#this script also generates a plot from the 100 selected Bxb1 helix sequences that best match the enriched patterns
#written by Jeff Miller for Sangamo Therapeutics

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


def format_tetrapeptide(tetrapeptide, pos1, pos2, pos3, pos4, peptide_length):
    output_string = ''
    for pos in range(peptide_length):
        if pos == pos1 or pos == pos2 or pos == pos3 or pos == pos4:
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

def pattern_create(seq1, seq2):
    pattern = ''
    for pos in range(len(seq1)):
        if seq1[pos] == seq2[pos]:
            pattern += seq1[pos]
        else:
            pattern += 'X'
    return pattern

if len(sys.argv) < 2:
    print('please give the input filename as a command line argument')
else:
    read_threshold = 2
    tetrapeptide_count_threshold = 3
    data_filename = sys.argv[1]
    truncated_filename = data_filename[:-4]  
    all_peptide_list = []
    above_threshold_peptide_list = []
    data_file = open(data_filename, 'r')
    for raw_line in data_file:
        line = raw_line.strip()
        if line:
            DNA, raw_DNA_count, peptide = line.split('\t')
            DNA_count = int(raw_DNA_count)
            all_peptide_list.append(peptide)
            if DNA_count >= read_threshold:
                above_threshold_peptide_list.append(peptide)

    all_peptide_dict = {}
    tetrapeptide_dict = {}
    total_peptides = len(all_peptide_list)
    total_above_threshold_peptides = len(above_threshold_peptide_list)
    peptide_length = len(all_peptide_list[0])
    exclude_pos_list = [4]
    peptide_offset = 231

    #intialize dictionaries
    for pos in range(peptide_length):
        all_peptide_dict[pos] = {}
        for residue in residue_list:
            all_peptide_dict[pos][residue] = 0
    for pos1 in range(peptide_length):
        for pos2 in range(peptide_length):
            for pos3 in range(peptide_length):
                for pos4 in range(peptide_length):
                    if pos2 > pos1 and pos3 > pos2 and pos4 > pos3 and pos1 not in exclude_pos_list and pos2 not in exclude_pos_list and pos3 not in exclude_pos_list and pos4 not in exclude_pos_list:
                        tetrapeptide_dict[(pos1,pos2,pos3,pos4)] = {}
                        for tetrapeptide in tetrapeptide_list:
                            tetrapeptide_dict[(pos1,pos2,pos3,pos4)][tetrapeptide] = 0




    for peptide in all_peptide_list:
        for pos in range(peptide_length):
            all_peptide_dict[pos][peptide[pos]] += 1
    for peptide in above_threshold_peptide_list:
        for pos1 in range(peptide_length):
            for pos2 in range(peptide_length):
                for pos3 in range(peptide_length):
                    for pos4 in range(peptide_length):
                        if pos2 > pos1 and pos3 > pos2 and pos4 > pos3 and pos1 not in exclude_pos_list and pos2 not in exclude_pos_list and pos3 not in exclude_pos_list and pos4 not in exclude_pos_list:
                            tetrapeptide = peptide[pos1] + peptide[pos2] + peptide[pos3] + peptide[pos4]
                            tetrapeptide_dict[(pos1,pos2,pos3,pos4)][tetrapeptide] += 1
    output_line = '%s\t%d\t' %(truncated_filename, total_above_threshold_peptides)

    residue_frequency_dict = {}
    for pos in range(peptide_length):
        residue_frequency_dict[pos] = {}
        for res in residue_list:
            residue_frequency_dict[pos][res] = float(all_peptide_dict[pos][res]) / float(total_peptides)
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
                                            formatted_tetrapeptide = format_tetrapeptide(tetrapeptide, pos1, pos2, pos3, pos4, peptide_length)
                                            observed_tetrapeptide = tetrapeptide_dict[(pos1, pos2, pos3, pos4)][tetrapeptide]
                                            expected_tetrapeptide = float(total_above_threshold_peptides) * (residue_frequency_dict[pos1][res1] * residue_frequency_dict[pos2][res2] * residue_frequency_dict[pos3][res3] * residue_frequency_dict[pos4][res4])
                                            expected_fraction = (residue_frequency_dict[pos1][res1] * residue_frequency_dict[pos2][res2] * residue_frequency_dict[pos3][res3] * residue_frequency_dict[pos4][res4])
                                            enrichment = float(observed_tetrapeptide) / expected_tetrapeptide
                                            if enrichment > 1:
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
    pattern_filename = truncated_filename + '_4res_patterns.txt'
    pattern_file = open(pattern_filename, 'w')
    filtered_pattern_list = sorted_enrichment_list[:100]
    final_pattern_list = []
    for item in filtered_pattern_list:
        pval = item[1][0]
        #print(pval)
        corrected_pval = pval_correction_dict[pval]
        #print(corrected_pval)
        if corrected_pval < 0.05 and item[1][2] >= tetrapeptide_count_threshold:
            #print('%s\t%6.2e\t%5.2f\t%d' %(item[0], corrected_pval, item[1][1], item[1][2]))
            pattern_file.write('%s\t%s\t%6.2e\t%5.2f\t%d\n' %(truncated_filename, item[0], corrected_pval, item[1][1], item[1][2]))
            final_pattern_list.append(item)

    data_file.close()
    peptide_pattern_score_dict = {}
    for peptide in above_threshold_peptide_list:
        pattern_match_score = 0
        for pattern in final_pattern_list:
            if pattern_match(peptide, pattern[0]) == 1:
                if pattern[1][0] > 0:
                    pattern_match_score += -1.0 * math.log(pattern[1][0], 10)
                else:
                    pattern_match_score += 312.0 # p values below abount 10^-312 are set to 0 and cause math.log to throw and error
        peptide_pattern_score_dict[peptide] = pattern_match_score

    sorted_peptide_score_list = sorted(peptide_pattern_score_dict.items(), key=lambda x:x[1], reverse=True)
    plot_peptide_list = []
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

    
    for item in sorted_peptide_score_list[:100]:
        #print('%s\t%5.2f' %(item[0], item[1]))
        if item[1] > 0:
            #output_file.write('%s\t%5.2f\n' %(item[0], item[1]))
            plot_peptide_list.append(item[0])
    if len(plot_peptide_list) > 0:
        for peptide in plot_peptide_list:
            for pos in range(peptide_length):
                residue_pos_dict[peptide[pos]][pos] += 1
        for pos in range(peptide_length):
            for res in residue_list:
                residue_pos_prob_dict[res][pos + peptide_offset] = float(residue_pos_dict[res][pos]) / float(len(plot_peptide_list))
        output_line += 'Information content at each position:\t'
        for pos in range(peptide_length):
            IC_value = IC(residue_pos_prob_dict, pos + peptide_offset)
            output_line += '%5.2f\t' %IC_value
            #logfile_output_line += '%5.2f\t' %IC_value
            for res in residue_list:
                IC_scaled_residue_pos_prob_dict[res][pos + peptide_offset] =  residue_pos_prob_dict[res][pos + peptide_offset] * IC_value   
        output_line += '%d\t' %len(plot_peptide_list)
        print(output_line)
        df = pd.DataFrame(IC_scaled_residue_pos_prob_dict)
        color_scheme = {
        'A' : [0, 0, 0.8],
        'C' : [0, 0, 0.8],
        'D' : [1, 0, 1],
        'E' : [1, 0, 1],
        'F' : [0, 0, 0],
        'G' : [1, .7, 0],
        'H' : [.5, .22, .05],
        'I' : [0, 0, 0.8],
        'K' : [0.8, 0, 0],
        'L' : [0, 0, 0.8],
        'M' : [0, 0, 0.8],
        'N' : [.6, 0, 0.5],
        'P' : [.9, .9, 0],
        'Q' : [.6, 0, 0.5],
        'R' : [0.8, 0, 0],
        'S' : [0, .5, 0],
        'T' : [0, .5, 0],
        'V' : [0, 0, 0.8],
        'W' : [0, 0, 0],
        'Y' : [0, 0, 0]}

        ww_logo = logomaker.Logo(df,
                             #color_scheme='weblogo_protein',
                             #color_scheme='NajafabadiEtAl2017',
                             color_scheme=color_scheme,
                             vpad= 0.01,
                             width=.9,
                             figsize=[8.0, 8.0])
        ww_logo.style_single_glyph(c='L', p=235, color=[0.85, 0.85, 0.85])
        ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
        #ww_logo.highlight_position(p=4, color='gold', alpha=.5)
        ww_logo.ax.set_ylabel('information (bits)')
        ww_logo.ax.text(peptide_offset, 4.5, '%s %d total peptides' %(truncated_filename, total_peptides), fontsize=12)
        ww_logo.ax.text(peptide_offset, -0.3, '%d plotted' %len(plot_peptide_list), fontsize=9)
                                                    
        ww_logo.ax.set_xlim([peptide_offset - 0.5, peptide_offset + peptide_length - 0.5])
        plot_filename = truncated_filename + '_pattern_matches.png'
        plt.savefig(plot_filename)
