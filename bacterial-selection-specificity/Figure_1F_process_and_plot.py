#this script analyzses the results of bacterial specificity selection
#input sequences have a 3bp region in the left halfsite randomized
#the reverse complement of this randomized region is in the right halfsite
#the left halfsite must match the reverse complement of the right halfsite for the sequence read to be counted
#output is a sequence motif plot that is the average of all replicates 
#written by Jeff Miller for Sangamo Therapeutics

import sys
from operator import itemgetter
import time
import logomaker 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math


BASES = ['A', 'C', 'G', 'T']
TRIPLETS = [] #list of DNA triplet sequences
for base1 in BASES:
    for base2 in BASES:
        for base3 in BASES:
            triplet = base1 + base2 + base3 
            TRIPLETS.append(triplet)

#caluclates information content of a given position based on the frequencies of each base
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
    print('please use a command line arguments with the degenerate template sequence and the list of .fq input filenames')
else:
    #read input file and convert into a list of filenames for input fastq files
    template_filename = sys.argv[1] #template sequence is oligo_NNN_library_3-6-2024.txt for this example
    template_file = open(template_filename, 'r')
    randomized_length = 3
    for raw_line in template_file:
        line = raw_line.strip()
        if line:
            template_sequence = line
            left_randomized_region_start = template_sequence.find('NNN')
            right_randomized_region_start = template_sequence[left_randomized_region_start + 1:].find('NNN') + left_randomized_region_start + 1
    input_filename_list_filename = sys.argv[2] 
    input_filename_list_file = open(input_filename_list_filename, 'r')
    input_filename_list = []
    base_pos_prob_dict = {} #dictionary of proportions of each base at each position in the sequence
    replicate = 0 
    for line in input_filename_list_file:
        line = line.strip()
        if line:
            print(line)
            input_filename_list.append(line)
            base_pos_prob_dict[replicate] = {}
            for base in BASES:
                base_pos_prob_dict[replicate][base] = {}
                for pos in range(3):
                    base_pos_prob_dict[replicate][base][pos] = 0.0
            replicate += 1

    helix = input_filename_list[10:17]
    replicate = 0
    #process the sequence file for each replicate
    for input_filename in input_filename_list:
        base_pos_dict = {}
        ct = 0
        total = 0
        length104 = 0
        length146 = 0
        constant_region_correct = 0
        randomized_region_inverted_repeat = 0
        constant_5p = template_sequence[:left_randomized_region_start]
        left_rand = template_sequence[left_randomized_region_start:left_randomized_region_start + randomized_length]
        constant_middle = template_sequence[left_randomized_region_start + randomized_length:right_randomized_region_start] 
        right_rand = template_sequence[right_randomized_region_start:right_randomized_region_start + randomized_length]
        constant_3p = template_sequence[right_randomized_region_start + randomized_length:]


        read_dict = {}
        input_file = open(input_filename, 'r')
        triplet_dict = {}
        for triplet in TRIPLETS:
            triplet_dict[triplet] = 0
        randomized_region_dict = {} #dictionary of each unique expected DNA sequence read and number of times that read was observed
        data_end = len(template_sequence) + 4

        #this loop reads in the fastq or fq file. Raw mode means 4 bp randomer needs to be removed
        for raw_line in input_file:
            line = raw_line.strip()
            if line:
                ct += 1
                if ct%4 == 1: #fastq files has four lines for each sequence read- first line is name of sequence read
                    seq_name = line
                if ct%4 == 2: #second line is actual sequence read (basecalls)
                    total += 1
                    seq = line #line[:data_end]
                    if len(seq) == 108: #4 bp randomer hasn't been trimmed yet
                        length104 += 1
                    if len(seq) == 150:
                        length146 += 1
                    if seq not in read_dict:
                        read_dict[seq] = 1
                    else:
                        read_dict[seq] += 1


        sorted_unique_read_list = sorted(read_dict.items(), key=lambda x:x[1], reverse=True) #sort sequence dictionary by read count with most reads first
        for base in BASES:
            base_pos_dict[base] = {}
            for pos in range(3):
                base_pos_dict[base][pos] = 0
        for item in sorted_unique_read_list:
            trimmed_read = item[0][4:]
            seq = trimmed_read
            if seq[:len(constant_5p)] == constant_5p and seq[left_randomized_region_start + randomized_length:right_randomized_region_start] == constant_middle and seq[-len(constant_3p):] == constant_3p: #only process files where constant region matches expectations
                constant_region_correct += item[1]
                left_randomized_region = seq[left_randomized_region_start:left_randomized_region_start + randomized_length]
                right_randomized_region = seq[right_randomized_region_start:right_randomized_region_start + randomized_length]
                #print('%s\t%s' %(left_randomized_region, right_randomized_region))
                if left_randomized_region == reverse_complement(right_randomized_region):
                    randomized_region_inverted_repeat += item[1]
                    triplet_dict[left_randomized_region] += item[1]
                    for pos in range(3):
                        base_pos_dict[left_randomized_region[pos]][pos] += item[1]

        for pos in range(3):
            for base in BASES:
                base_pos_prob_dict[replicate][base][pos] = float(base_pos_dict[base][pos]) / float(randomized_region_inverted_repeat)
        replicate += 1



    #averages the position weight matrices for all replicates
    avearage_base_pos_prob_dict = {}
    for base in BASES:
        avearage_base_pos_prob_dict[base] = {}
    for pos in range(3):
        for base in BASES:
            cumulative_base_pos_prob_dict = 0.0
            for replicate in range(len(input_filename_list)):
                cumulative_base_pos_prob_dict += base_pos_prob_dict[replicate][base][pos]
            avearage_base_pos_prob_dict[base][pos] = cumulative_base_pos_prob_dict / float(len(input_filename_list))

    IC_scaled_base_pos_prob_dict = {}
    for base in BASES:
        IC_scaled_base_pos_prob_dict[base] = {}
        for pos in range(3):
            IC_scaled_base_pos_prob_dict[base][pos] = 0.0

    for pos in range(3):
        IC_value = IC(avearage_base_pos_prob_dict, pos)
        for base in BASES:
            IC_scaled_base_pos_prob_dict[base][pos - 11] =  avearage_base_pos_prob_dict[base][pos] * IC_value   


    df = pd.DataFrame(IC_scaled_base_pos_prob_dict)

    ww_logo = logomaker.Logo(df,
                         color_scheme='classic',
                         vpad= 0.01,
                         width=.9,
                         figsize=[3.0, 3.0])
    ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
    ww_logo.ax.set_ylabel('Information (bits)')
    ww_logo.ax.set_yticks([0, 1, 2])
    ww_logo.ax.set_xlim([- 11.5, -8.5])
    ww_logo.ax.set_ylim([0, 2])
    plot_filename = input_filename_list_filename[:-13] +  '.pdf' 
    plt.savefig(plot_filename)

                    
