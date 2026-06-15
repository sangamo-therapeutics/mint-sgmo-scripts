# ExtD_Figure_2F_process_and_plot.py
# this script analyzses the results of bacterial specificity selection
# input sequences have a 2bp region corresponding to a Bxb1 loop target in the left halfsite randomized
# the reverse complement of this randomized region is in the right halfsite
# the left halfsite must match the reverse complement of the right halfsite for the sequence read to be counted
# output is a sequence motif plot that is the average of all replicates
# written by Jeff Miller for Sangamo Therapeutics

"""This script can be invoked with the command

    $  python Figure_1F_process_and_plot.py /path/to/oligo_NN_library_loop.txt \
            /path/to/ExtD_FIGURE_2F_YRGSLP_fq_files.txt \
            /path/to/location/of/fq_files

    The third argument is optional, but if not supplied, data files must be in the same directory as the
    ExtD_FIGURE_2F_YRGSLP_fq_files.txt

    The output .pdf containing the figure will be written to a subdirectory "figure_output"

"""

import math
import sys
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import pandas as pd

BASES = ['A', 'C', 'G', 'T']
DIMERS = []
for base1 in BASES:
    for base2 in BASES:
        dimer = base1 + base2
        DIMERS.append(dimer)


# caluclates information content of a given position based on the frequencies of each base
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
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', ' ': ' ', '\n': '\n', \
                  'B': 'V', 'D': 'H', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', \
                  'R': 'Y', 'S': 'S', 'V': 'B', 'W': 'W', 'Y': 'R'}  # dictionary of WC basepairs
    seq = list(forward_sequence)  # turns sequence into list
    seq.reverse()
    rev_comp = []
    for x in seq:
        rev_comp.append(complement[x])
    reverse_sequence = ''
    for x in rev_comp:  # turns processed list back into string
        reverse_sequence = reverse_sequence + x
    return (reverse_sequence)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError('please use a command line arguments with the degenerate template '
                         'sequence and the list of .fq input filenames')

    # read input file and convert into a list of filenames for input fastq files
    template_filename = sys.argv[1]  # template sequence is oligo_NNN_library_3-6-2024.txt for this example
    print(f"---Template file is {template_filename}")
    # template_filename = r'C:\Users\jeshleman\PycharmProjects\mint-sgmo-scripts\data\figure_files\oligo_NN_library_loop.txt'
    input_filename_list_filename = sys.argv[2]
    print(f"-----Data source file list is {input_filename_list_filename}")
    input_filename_list_filename = r'C:\Users\jeshleman\PycharmProjects\mint-sgmo-scripts\scripts\figure_maker_scripts\ExtD_FIGURE_2F_YRGSLP_fq_files.txt'
    try:
        data_location = sys.argv[3]
        print(f"------We will look for data in {data_location}")
    except IndexError:
        data_location = Path(input_filename_list_filename).parent
        print(f"------We will use the default data location {data_location}")
    print(f"-----------Data location is {data_location}")
    data_location = Path(data_location)

    template_file = Path(template_filename)
    randomized_length = 2

    template_sequence = template_file.read_text().strip()

    if not template_sequence:
        raise ValueError(f"Template file {template_file.name} did not contain a sequence")

    left_randomized_region_start = template_sequence.find('NN')
    right_randomized_region_start = template_sequence[left_randomized_region_start + 1:].find(
        'NN') + left_randomized_region_start + 1

    input_filename_list_file = open(input_filename_list_filename, 'r')
    input_filename_list = []
    base_pos_prob_dict = {}  # dictionary of proportions of each base at each position in the sequence
    replicate = 0

    for line in input_filename_list_file:
        line = line.strip()
        if line:
            print(line)
            input_filename_list.append(line)
            base_pos_prob_dict[replicate] = {}
            for base in BASES:
                base_pos_prob_dict[replicate][base] = {}
                for pos in range(2):
                    base_pos_prob_dict[replicate][base][pos] = 0.0
            replicate += 1
    print()
    loop_seq = input_filename_list[0][
               15:21]  # the amino acid sequence of the loop variant is found within the filename

    # dictionaries to store the data for each replicate
    base_pos_dict = {}
    base_pos_prob_dict = {}
    IC_scaled_base_pos_prob_dict = {}

    replicate = 0
    # each filename is considered a separate replicate and processed by the loop below
    for input_filename in input_filename_list:
        dimer_dict_total = {}
        dimer_dict_recombined = {}
        dimer_dict_recombined_inverted_repeat = {}
        for dimer in DIMERS:
            dimer_dict_total[dimer] = 0
            dimer_dict_recombined[dimer] = 0
            dimer_dict_recombined_inverted_repeat[dimer] = 0

        base_pos_dict = {}
        ct = 0
        total = 0
        constant_region_correct = 0
        randomized_region_inverted_repeat = 0
        constant_5p = template_sequence[:left_randomized_region_start]
        left_rand = template_sequence[left_randomized_region_start:left_randomized_region_start + randomized_length]
        constant_middle = template_sequence[
                          left_randomized_region_start + randomized_length:right_randomized_region_start]
        right_rand = template_sequence[right_randomized_region_start:right_randomized_region_start + randomized_length]
        constant_3p = template_sequence[right_randomized_region_start + randomized_length:]

        read_dict = {}
        # input_file = open(input_filename, 'r')
        input_file = open(data_location / input_filename, 'r')
        dimer_dict = {}
        for dimer in DIMERS:
            dimer_dict[dimer] = 0
        randomized_region_dict = {}  # dictionary of each unique expected DNA sequence read and number of times that read was observed
        data_end = len(template_sequence) + 4

        # this loop reads in the fastq or fq file. Raw mode means 4 bp randomer needs to be removed
        for raw_line in input_file:
            line = raw_line.strip()
            if line:
                ct += 1
                if ct % 4 == 1:  # fastq files has four lines for each sequence read- first line is name of sequence read
                    seq_name = line
                if ct % 4 == 2:  # second line is actual sequence read (basecalls)
                    total += 1
                    seq = line  # line[:data_end]
                    if seq not in read_dict:
                        read_dict[seq] = 1
                    else:
                        read_dict[seq] += 1

        sorted_unique_read_list = sorted(read_dict.items(), key=lambda x: x[1],
                                         reverse=True)  # sort sequence dictionary by read count with most reads first
        for item in sorted_unique_read_list:
            trimmed_read = item[0][:]
            seq = trimmed_read
            if seq[
               :len(
                   constant_5p)] == constant_5p:  # and seq[left_randomized_region_start + randomized_length:right_randomized_region_start] == constant_middle and seq[-len(constant_3p):] == constant_3p: #only process files where constant region matches expectations
                constant_region_correct += item[1]  # 1 retains duplicate reads
                left_randomized_region = seq[
                                         left_randomized_region_start:left_randomized_region_start + randomized_length]
                right_randomized_region = seq[
                                          right_randomized_region_start:right_randomized_region_start + randomized_length]
                dimer_dict_total[left_randomized_region] += item[1]  # retains duplicate reads
                # all constant regions of the sequence read must match the expected sequence or the read is filtered out
                if seq[:len(constant_5p)] == constant_5p and seq[
                                                             left_randomized_region_start + randomized_length:right_randomized_region_start] == constant_middle and seq[
                                                                                                                                                                    -len(
                                                                                                                                                                        constant_3p):] == constant_3p:
                    dimer_dict_recombined[left_randomized_region] += item[1]  # retains duplicate reads
                    # checks to make sure variable sequence in right half site is reverse complement of variable sequence in left halfsite
                    if left_randomized_region == reverse_complement(right_randomized_region):
                        randomized_region_inverted_repeat += item[1]
                        dimer_dict_recombined_inverted_repeat[left_randomized_region] += item[1]

        # initilize dictionaries
        for base in BASES:
            base_pos_dict[base] = {}
            for pos in range(2):
                base_pos_dict[base][pos] = 0
        base_pos_dict[replicate] = {}
        base_pos_prob_dict[replicate] = {}
        IC_scaled_base_pos_prob_dict[replicate] = {}
        for base in BASES:
            base_pos_dict[replicate][base] = {}
            base_pos_prob_dict[replicate][base] = {}
            IC_scaled_base_pos_prob_dict[replicate][base] = {}
            for pos in range(2):
                base_pos_dict[replicate][base][pos] = 0
                base_pos_prob_dict[replicate][base][pos] = 0.0
                IC_scaled_base_pos_prob_dict[replicate][base][pos] = 0.0

        output_line = '%s\t%d\t%d\t%d\t%5.2f\t%d\t' % (
            input_filename, total, len(sorted_unique_read_list), constant_region_correct,
            100.0 * float(constant_region_correct) / float(total), randomized_region_inverted_repeat)
        for dimer in DIMERS:
            if dimer_dict_total[dimer] > 0:
                output_line += '%5.4f\t' % (float(dimer_dict_recombined_inverted_repeat[dimer]) / float(
                    randomized_region_inverted_repeat))
            else:
                output_line += 'NA\t'
            for pos in range(2):
                base_pos_dict[replicate][dimer[pos]][pos] += dimer_dict_recombined_inverted_repeat[dimer]
        print(output_line)
        for pos in range(2):
            for base in BASES:
                base_pos_prob_dict[replicate][base][pos] = float(base_pos_dict[replicate][base][pos]) / float(
                    randomized_region_inverted_repeat)
        for pos in range(2):
            for base in BASES:
                output_line += '%5.2f\t' % (base_pos_prob_dict[replicate][base][pos])
        output_line = 'Information content at each position:\t'
        for pos in range(randomized_length):
            IC_value = IC(base_pos_prob_dict[replicate], pos)
            output_line += '%5.2f\t' % IC_value
            for base in BASES:
                IC_scaled_base_pos_prob_dict[replicate][base][pos] = base_pos_prob_dict[replicate][base][pos] * IC_value
        print(output_line)
        replicate += 1

    # averages the position weight matrices for all replicates
    avearage_base_pos_prob_dict = {}
    for base in BASES:
        avearage_base_pos_prob_dict[base] = {}
    for pos in range(2):
        for base in BASES:
            cumulative_base_pos_prob_dict = 0.0
            for replicate in range(len(input_filename_list)):
                cumulative_base_pos_prob_dict += base_pos_prob_dict[replicate][base][pos]
            avearage_base_pos_prob_dict[base][pos] = cumulative_base_pos_prob_dict / float(len(input_filename_list))

    IC_scaled_base_pos_prob_dict = {}
    for base in BASES:
        IC_scaled_base_pos_prob_dict[base] = {}
        for pos in range(2):
            IC_scaled_base_pos_prob_dict[base][pos] = 0.0

    for pos in range(2):
        IC_value = IC(avearage_base_pos_prob_dict, pos)  # calculated information content for motif plot
        for base in BASES:
            IC_scaled_base_pos_prob_dict[base][pos] = avearage_base_pos_prob_dict[base][pos] * IC_value

    df = pd.DataFrame(IC_scaled_base_pos_prob_dict)

    ww_logo = logomaker.Logo(df,
                             color_scheme='classic',
                             vpad=0.01,
                             width=.9,
                             figsize=[3.0, 3.0])
    ww_logo.style_xticks(anchor=0, spacing=1, rotation=0)
    ww_logo.ax.set_ylabel('information (bits)')
    ww_logo.ax.set_yticks([0, 1.0, 2])
    ww_logo.ax.text(-0.2, 2.1, '%s with %d sequences' % (loop_seq, randomized_region_inverted_repeat), fontsize=9)
    ww_logo.ax.set_xlim([- 0.5, 2 - 0.5])
    ww_logo.ax.set_ylim([0, 2])
    plot_filename = loop_seq + '.pdf'
    output_loc = data_location / "figure_output"
    output_loc.mkdir(exist_ok=True)
    plt.savefig(output_loc / plot_filename)
