#ExtD_Figure_2G_process_and_plot.py
#this script analyzses the results of a bxb1 specficity experiment with 6 randomized bases in each halfsite
#input is the expected sequence read including degenerate DNA bases at the randomized position and a list of fastq filenames with the actual reads for each replicate
#output is a plot of the motif of all 6 bp sequences (hexamers) that occur in all replicates
#written by Jeff Miller for Sangamo Therapeutics

import math
from operator import itemgetter
import logomaker 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import statistics
from operator import itemgetter

BASES = ['A', 'C', 'G', 'T']
HEXAMERS = [] #list of all possible 6 bp DNA sequences
for base1 in BASES:
    for base2 in BASES:
        for base3 in BASES:
            for base4 in BASES:
                for base5 in BASES:
                    for base6 in BASES:
                        hexamer = base1 + base2 + base3 + base4 + base5 + base6
                        HEXAMERS.append(hexamer)

#calculates information content at each position of a list of DNA sequences
def IC(base_frequency_dict, pos):
    cumulative_value = 0
    for base in BASES:
        base_frequency = base_frequency_dict[base][pos]
        if base_frequency != 0.0:
            cumulative_value -= base_frequency * math.log(base_frequency, 2)
    IC_value = 2 - cumulative_value
    return IC_value

# This function returns the reverse-complement of a DNA sequence 
def reverse_complement(forward_sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', ' ': ' ', '\n': '\n',\
                        'B': 'V', 'D': 'H', 'H': 'D', 'K': 'M', 'M': 'K', 'N':'N',\
                        'R': 'Y', 'S': 'S', 'V': 'B', 'W': 'W', 'Y': 'R'}       # dictionary of WC basepairs
        seq = list(forward_sequence)                            # turns sequence into list
        seq.reverse()
        rev_comp = []
        for x in seq:
                rev_comp.append(complement[x])
        reverse_sequence = ''
        for x in rev_comp:                                      # turns processed list back into string
                reverse_sequence = reverse_sequence + x
        return(reverse_sequence)

if len(sys.argv) < 3:
    print('please use a command line arguments with the degenerate template sequence and a list of input file names')
else:
    ct = 0
    total = 0
    seq_length = 6
    unique_read_dict = {}
    constant_region_correct = 0
    randomized_region_inverted_repeat = 0

    #read expected sequence template with N at each randomized position in the library of partially randomized DNA sequences
    template_filename = sys.argv[1]
    template_file = open(template_filename, 'r')
    for raw_line in template_file:
        line = raw_line.strip()
        if line:
            template_sequence = line
            if template_sequence.find('NNNNNN') != -1:
                randomized_length = 6
                left_randomized_region_start = template_sequence.find('NNNNNN')
                right_randomized_region_start = template_sequence[left_randomized_region_start + 1:].find('NNNNNN') + left_randomized_region_start + 1
    constant_5p = template_sequence[:left_randomized_region_start]
    left_rand = template_sequence[left_randomized_region_start:left_randomized_region_start + randomized_length]
    constant_middle = template_sequence[left_randomized_region_start + randomized_length:right_randomized_region_start] 
    right_rand = template_sequence[right_randomized_region_start:right_randomized_region_start + randomized_length]
    constant_3p = template_sequence[right_randomized_region_start + randomized_length:]

    #read second input file and convert into a list of filenames for input fastq files
    input_filename_list_filename = sys.argv[2] 
    input_filename_list_file = open(input_filename_list_filename, 'r')
    input_filename_list = []
    base_pos_prob_dict = {}
    replicate = 0
    fraction_list_dict = {}
    average_fraction_dict = {}
    for hexamer in HEXAMERS:
        fraction_list_dict[hexamer] = []
        average_fraction_dict[hexamer] = 0.0
    for line in input_filename_list_file:
        line = line.strip()
        if line:
            print(line)
            input_filename_list.append(line)
    hairpin_seq = input_filename_list[0][15:27] #the amino acid sequence of the loop variant is found within the filename 

    for input_filename in input_filename_list: 
        input_file = open(input_filename, 'r')
        randomized_region_dict = {} #dictionary of each unique expected DNA sequence read and number of times that read was observed
        data_end = len(template_sequence) + 4
        #this loop reads in the fastq or fq file. 
        for raw_line in input_file:
            line = raw_line.strip()
            if line:
                ct += 1
                if ct%4 == 1: #fastq files has four lines for each sequence read- first line is name of sequence read
                    seq_name = line
                if ct%4 == 2: #second line is actual sequence read (basecalls)
                    total += 1
                    seq = line[4:data_end] #remove randomized 4 bp at start of read
                    if seq[:len(constant_5p)] == constant_5p and seq[left_randomized_region_start + randomized_length:right_randomized_region_start] == constant_middle and seq[-len(constant_3p):] == constant_3p: #only process files where constant region matches expectations
                        constant_region_correct += 1
                        left_randomized_region = seq[left_randomized_region_start:left_randomized_region_start + randomized_length]
                        right_randomized_region = seq[right_randomized_region_start:right_randomized_region_start + randomized_length]
                        if left_randomized_region == reverse_complement(right_randomized_region):
                            randomized_region_inverted_repeat += 1
                            if left_randomized_region not in randomized_region_dict: #randomized regions that matches intended randomization added to dictionary
                                randomized_region_dict[left_randomized_region] = 1
                            else:   
                                randomized_region_dict[left_randomized_region] += 1
        sorted_randomized_list = sorted(randomized_region_dict.items(), key=lambda x:x[1], reverse=True) #sort sequence dictionary by read count with most reads first
        for item in sorted_randomized_list:
            fraction_list_dict[item[0]].append(float(item[1])/float(randomized_region_inverted_repeat))

    filtered_hexamer_list = []
    for hexamer in HEXAMERS:
        if len(fraction_list_dict[hexamer]) == len(input_filename_list): #only analyzes sequences that occur in all replicates
            average_fraction_dict[hexamer] = statistics.mean(fraction_list_dict[hexamer])
            filtered_hexamer_list.append(hexamer)

    #initialize dictionaries
    base_pos_dict = {} #dictionary of counts of each base at each position within the hexamer
    base_pos_prob_dict = {} #normalized version of previous dictionary
    IC_scaled_base_pos_prob_dict = {}
    for base in BASES:
        base_pos_dict[base] = {}
        base_pos_prob_dict[base] = {}
        IC_scaled_base_pos_prob_dict[base] = {}
        for pos in range(seq_length):
            base_pos_dict[base][pos] = 0
            base_pos_prob_dict[base][pos] = 0.0
            IC_scaled_base_pos_prob_dict[base][pos] = 0.0

    cumulative_motif_fraction = 0.0 #this corrects for sequences that passed the filtering step
    for seq in filtered_hexamer_list:
        cumulative_motif_fraction += average_fraction_dict[seq]
        for pos in range(seq_length):
            base_pos_dict[seq[pos]][pos] += average_fraction_dict[seq]
    for pos in range(seq_length):
        for base in BASES:
            base_pos_prob_dict[base][pos] = base_pos_dict[base][pos] / cumulative_motif_fraction #applies correction
    for pos in range(seq_length):
        IC_value = IC(base_pos_prob_dict, pos)
        for base in BASES:
            IC_scaled_base_pos_prob_dict[base][pos - 18] =  base_pos_prob_dict[base][pos] * IC_value   

    #generate DNA sequence logo plot
    plot_filename = hairpin_seq + '_average.pdf'
    df = pd.DataFrame(IC_scaled_base_pos_prob_dict)
    ww_logo = logomaker.Logo(df,
         color_scheme='classic',
         vpad= 0.01,
         width=.9,
         figsize=[8.0, 8.0])
    ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
    ww_logo.ax.set_ylabel('information (bits)')
    ww_logo.ax.text(-18, 2.1, '%s ' %(hairpin_seq), fontsize=12)
    ww_logo.ax.set_xlim([- 18.5, -12.5])
    ww_logo.ax.set_ylim([0, 2])
    plt.savefig(plot_filename)


                        
