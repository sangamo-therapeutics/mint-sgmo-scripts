#this script identifies motifs in a list of enriched patterns from a bxb1 loop selection 
#input is a sorted list of enriched 4 residue patterns
#output is a plot of the sequences that best match the enriched patters in each motif as well as details about each motif
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


def pattern_match(peptide, pattern):
    match = 1
    for pos in range(len(pattern)):
        if pattern[pos] != 'X' and peptide[pos] != 'X' and pattern[pos] != peptide[pos]:
            match = 0
    return match

def pattern_compare(pattern1, pattern2):
    matches = 0
    mismatches = 0
    for pos in range(len(pattern1)):
        if pattern1[pos] != 'X' and pattern2[pos] != 'X' and pos != 8:
            if pattern1[pos] == pattern2[pos]:
                matches += 1
            else:
                mismatches += 1
    return (matches, mismatches)

def IC(residue_frequency_dict, pos):
    cumulative_value = 0
    for residue in residue_list:
        residue_frequency = residue_frequency_dict[residue][pos]
        if residue_frequency != 0.0:
            cumulative_value -= residue_frequency * math.log(residue_frequency, 2)
    IC_value = 4.32193 - cumulative_value
    return IC_value



if len(sys.argv) < 2:
    print('please give the input filename as a command line argument')
else:
    peptide_length = 6
    peptide_offset = 154
    peptide_list = []
    pattern_list = []
    data_filename = sys.argv[1]
    truncated_filename = data_filename[:-18]
    peptide_filename = truncated_filename + '.txt'
    peptide_file = open(peptide_filename, 'r')
    for raw_line in peptide_file:
        line = raw_line.strip()
        DNA, raw_count, peptide = line.split()
        peptide_list.append(peptide)
    print(len(peptide_list))
    print(truncated_filename)
    data_file = open(data_filename, 'r')
    for raw_line in data_file:
        line = raw_line.strip()
        if line and line[0] != '#':
            pattern_filename, pattern, raw_pval, raw_enrich, raw_count = line.split()
            pval = float(raw_pval)
            pattern_list.append((pattern, pval))
    sorted_pattern_list = sorted(pattern_list, key=lambda x:x[1], reverse=False)

    in_motif_list = []
    not_in_motif_list = []
    for pos in range(len(sorted_pattern_list)):
        not_in_motif_list.append(pos)
    motif_count = 0
    motif_dict = {}
    motif_size_at_start_of_round = 0
    motif_seed = 0 #strongest pattern
    rounds = 0
    while motif_seed != not_in_motif_list[-1]:
        motif_dict[motif_count] = [motif_seed]
        in_motif_list.append(motif_seed)
        not_in_motif_list.remove(motif_seed)
        #print(motif_seed)
        while len(motif_dict[motif_count]) != motif_size_at_start_of_round:
            rounds += 1
            motif_size_at_start_of_round = len(motif_dict[motif_count])
            unmatched_pattern_positions = []
            for pos in range(len(sorted_pattern_list)):
                best_matches = 0
                if pos not in in_motif_list:
                    for item in motif_dict[motif_count]:
                        pattern1 = sorted_pattern_list[item][0]
                        pattern2 = sorted_pattern_list[pos][0]
                        matches = pattern_compare(pattern1, pattern2)[0]
                        if matches > best_matches:
                            best_matches = matches
                    pattern1 = sorted_pattern_list[motif_seed][0]
                    pattern2 = sorted_pattern_list[pos][0]
                    seed_mismatches = pattern_compare(pattern1, pattern2)[1]
                if best_matches >= 3 and seed_mismatches <= 1:
                    motif_dict[motif_count].append(pos)
                    in_motif_list.append(pos)
                    not_in_motif_list.remove(pos)
                elif pos not in in_motif_list:
                    unmatched_pattern_positions.append(pos)
        motif_count += 1
        if len(not_in_motif_list) == 0:
            break
        motif_seed = not_in_motif_list[0]
        motif_dict[motif_count] = [motif_seed]

    #print()
    motif_pos = 0
    output_motif_num = 1
    while motif_pos < motif_count:
        if len(motif_dict[motif_pos]) >= 5:
            motif_filename = truncated_filename + '_motif%d.txt' %output_motif_num
            motif_file = open(motif_filename, 'w')
            for item in motif_dict[motif_pos]:
                print('%d\t%s\t%5.2e' %(item, sorted_pattern_list[item][0], sorted_pattern_list[item][1]))
                motif_file.write('%d\t%s\t%5.2e\n' %(item, sorted_pattern_list[item][0], sorted_pattern_list[item][1]))
            print()
            plot_filename = truncated_filename + '_motif%d.png' %output_motif_num
            
            output_motif_num += 1
            motif_peptides = []
            for peptide in peptide_list:
                include_in_motif_list = 0
                for item in motif_dict[motif_pos]:
                    if pattern_match(peptide, sorted_pattern_list[item][0]) == 1:
                        include_in_motif_list = 1
                        #print('%s\t%s' %(peptide, sorted_pattern_list[item][0]))
                if include_in_motif_list == 1:
                    motif_peptides.append(peptide)
                    motif_file.write('%s\n' %peptide)
            print(len(motif_peptides))
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

            for peptide in motif_peptides:
                for pos in range(peptide_length):
                    residue_pos_dict[peptide[pos]][pos] += 1
            for pos in range(peptide_length):
                for res in residue_list:
                    residue_pos_prob_dict[res][pos + peptide_offset] = float(residue_pos_dict[res][pos]) / float(len(motif_peptides))
            output_line = 'Information content at each position:\t'
            for pos in range(peptide_length):
                IC_value = IC(residue_pos_prob_dict, pos + peptide_offset)
                output_line += '%5.2f\t' %IC_value
                #logfile_output_line += '%5.2f\t' %IC_value
                for res in residue_list:
                    IC_scaled_residue_pos_prob_dict[res][pos + peptide_offset] =  residue_pos_prob_dict[res][pos + peptide_offset] * IC_value   
            #print(output_line)
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
            ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
            #ww_logo.highlight_position(p=4, color='gold', alpha=.5)
            ww_logo.ax.set_ylabel('information (bits)')
            ww_logo.ax.text(peptide_offset, 4.5, '%s - motif %d' %(truncated_filename, output_motif_num - 1), fontsize=12)
            ww_logo.ax.text(peptide_offset, -0.3, '%d peptides total: %d plotted' %(len(peptide_list), len(motif_peptides)), fontsize=9)
                                                            
            ww_logo.ax.set_xlim([peptide_offset - 0.5, peptide_offset + peptide_length - 0.5])
            
            plt.savefig(plot_filename)
            motif_file.close()
            
        motif_pos += 1
        
