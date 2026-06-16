#this script analyzses the results of a directed evolution experiment with a randomized Bxb1 helix
#input is the expected sequence read including degenerate DNA bases at the randomized position and an R1.fq fastq filename with the actual reads
#output file contains a line for each expected DNA sequence that encodes an open reading frame
#output is DNA sequence of randomized region, sequence count, translated peptide, number of peptides with a single mismatch, and peptides with a single amino acid mismatch and more than one DNA mismatch 
#written by Jeff Miller for Sangamo Therapeutics

import sys
from operator import itemgetter
import math

residue_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
BASES = ['A', 'C', 'G', 'T']
number_of_peptides_to_plot = 30

#translate DNA sequence into peptide sequence. Translates when DNA degeneracy doesn't make amino acid uncertain. Stop codon indicated by period.
def translate(DNA_seq): 
        genetic_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', \
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', \
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '.', 'TAG': '.', \
                'TGT': 'C', 'TGC': 'C', 'TGA': '.', 'TGG': 'W', \
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', \
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', \
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', \
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', \
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', \
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', \
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', \
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', \
                'AAN': 'X', 'ACN': 'T', 'AGN': 'X', 'ATN': 'X', \
                'CAN': 'X', 'CCN': 'P', 'CGN': 'R', 'CTN': 'L', \
                'GAN': 'X', 'GCN': 'A', 'GGN': 'G', 'GTN': 'V', \
                'TAN': 'X', 'TCN': 'S', 'TGN': 'X', 'TTN': 'X', \
                'ANA': 'X', 'ANC': 'X', 'ANG': 'X', 'ANT': 'X', \
                'GNA': 'X', 'CNC': 'X', 'CNG': 'X', 'CNT': 'X', \
                'TNA': 'X', 'GNC': 'X', 'GNG': 'X', 'GNT': 'X', \
                'TNA': 'X', 'TNC': 'X', 'TNG': 'X', 'TNT': 'X', \
                'NAA': 'X', 'NAC': 'X', 'NAG': 'X', 'NAT': 'X', \
                'NCA': 'X', 'NCC': 'X', 'NCG': 'X', 'NCT': 'X', \
                'NGA': 'X', 'NGC': 'X', 'NGG': 'X', 'NGT': 'X', \
                'NTA': 'X', 'NTC': 'X', 'NTG': 'X', 'NTT': 'X', \
                'NNA': 'X', 'NNC': 'X', 'NNG': 'X', 'NNT': 'X', \
                'NAN': 'X', 'NCN': 'X', 'NGN': 'X', 'NTN': 'X', \
                'ANN': 'X', 'CNN': 'X', 'GNN': 'X', 'TNN': 'X', \
                'NNN': 'X'}
        
        codon_position = 0
        peptide = ''
        while 1:
                codon = DNA_seq[codon_position * 3: codon_position * 3 + 3]
                if len(codon) < 3:break
                peptide = peptide + genetic_code[codon]
                codon_position = codon_position + 1
        return(peptide)

def IC(residue_frequency_dict, pos):
    cumulative_value = 0
    for residue in residue_list:
        residue_frequency = residue_frequency_dict[residue][pos]
        if residue_frequency != 0.0:
            cumulative_value -= residue_frequency * math.log(residue_frequency, 2)
    IC_value = 4.32193 - cumulative_value
    return IC_value



#finds position of first and last degenerate base in expected sequence. Returns tuple with first and last position of degenerate bases
def parse_degen_DNA_seq(seq):
    #print(seq)
    BASES = ['A', 'C', 'G', 'T']
    first_degen_base = -1
    last_degen_base = 1000
    for pos in range(len(seq)):
        if first_degen_base == -1 and seq[pos] not in BASES:
            first_degen_base = pos
        if last_degen_base == 1000 and seq[len(seq) - pos - 1] not in BASES:
            last_degen_base = len(seq) - pos - 1
    return (first_degen_base, last_degen_base)

def find_fixed_codons(seq):
    fixed_position_list = []
    for ct in range(int(len(seq) / 3)):
        if seq[(ct * 3)] in ['A', 'C', 'G', 'T'] and seq[(ct * 3 + 1)] in ['A', 'C', 'G', 'T'] and seq[(ct * 3 + 2)] in ['A', 'C', 'G', 'T']:
            fixed_position_list.append(ct)
    return fixed_position_list

#checks is sequence read is consistent with input randomized sequence- returns 1 for consistent and 0 for not consistent
def check_randomized_region(randomized_template, randomized_read):
    degen_dict = {'A':['A'], 'C':['C'], 'G':['G'], 'T':['T'],\
                  'N':['A', 'C', 'G', 'T'], 'K':['G', 'T'], 'S':['G', 'C']}
    correct_read = 1
    for pos in range(len(randomized_template)):
        if randomized_read[pos] not in degen_dict[randomized_template[pos]]:
            correct_read = 0
    return correct_read

#counts mismatches between two non-degenerate sequences of equal length- breaks after mismatch_threshold reached and returns mismatch_threshold
def count_mismatches(seq1, seq2, mismatch_threshold):
    mismatches = 0
    for pos in range(len(seq1)):
        if seq1[pos] != seq2[pos]:
            mismatches += 1
        if mismatches >= mismatch_threshold:
            break
    return mismatches


if len(sys.argv) < 3:
    print('please use a command line arguments with the degenerate template sequence, read count threshold, and the .fq input file name')
else:
    #read input file and convert into a list of filenames for input fastq files
    ct = 0
    total = 0
    read_count_threshold = 1 #set this higher for really noisy data
    constant_region_correct = 0
    template_filename = sys.argv[1]
    template_file = open(template_filename, 'r')
    for raw_line in template_file:
        line = raw_line.strip()
        if line:
            template_sequence = line
            randomized_start, randomized_end = parse_degen_DNA_seq(template_sequence)
    constant_5p = template_sequence[:randomized_start]
    constant_3p = template_sequence[randomized_end + 1:]
    randomized_template_region = template_sequence[randomized_start:randomized_end + 1]
    fixed_position_list = find_fixed_codons(randomized_template_region)
    randomized_size = randomized_end - randomized_start + 1
    peptide_length = int(randomized_size / 3)
    peptide_offset = 1
    input_filename = sys.argv[2] 
    raw_mode = 1
    if input_filename[-6:] == '_R1.fq':
        raw_mode = 1
        truncated_filename = input_filename[:-6]
    elif input_filename[-9:] == '_R1.fastq':
        raw_mode = 1
        truncated_filename = input_filename[:-9]
    else:
        truncated_filename = input_filename[:-3]
##    if raw_mode == 1:
##        print('raw mode- will trim first 4 bp')

        
    output_all_filename = truncated_filename + '_peptides_all.txt' #revmoes _R1.fq
    input_file = open(input_filename, 'r')
    output_all_file = open(output_all_filename, 'w')
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
                if raw_mode == 1:
                    seq = line[4:data_end]
                else:
                    seq = line
                if seq[:len(constant_5p)] == constant_5p and seq[-len(constant_3p):] == constant_3p: #only process files where constant region matches expectations
                    randomized_region = seq[len(constant_5p):len(constant_5p) + randomized_size]
                    if check_randomized_region(randomized_template_region, randomized_region) == 1: #checks randomized region to make sure it matches intended randomization scheme
                        constant_region_correct += 1
                        if randomized_region not in randomized_region_dict: #randomized regions that matches intended randomization added to dictionary
                            randomized_region_dict[randomized_region] = 1
                        else:   
                            randomized_region_dict[randomized_region] += 1

    print('%d total reads with %d matching expected pattern' %(total, constant_region_correct))
    sorted_randomized_list = sorted(randomized_region_dict.items(), key=lambda x:x[1], reverse=True) #sort sequence dictionary by read count with most reads first
    print('%d unique reads' %len(sorted_randomized_list))
    above_threshold_read_list = []
    for item in sorted_randomized_list:
        if item[1] >= read_count_threshold:
            above_threshold_read_list.append(item)
    print('%d unique reads above read count threshold' %len(above_threshold_read_list))
    filtered_read_list = [] 
    #rare reads that are very similar to extremely common reads are likely to be sequencing artifacts and are filtered out.
    #Filtering is accomplished by setting the entry in unique_dict for a given sequence to 0 if it is to be filtered out
    unique_dict = {}
    for pos in range(len(above_threshold_read_list)):
        unique_dict[pos] = 1
        seq = above_threshold_read_list[pos][0]
##    filter_start = time.time()
    for pos1 in range(len(above_threshold_read_list)):
        if above_threshold_read_list[pos1][1] > 50:
            for pos2 in range(len(above_threshold_read_list)):
                if pos2 > pos1:# and unique_dict[pos2] == 1:
                    read_count_ratio = above_threshold_read_list[pos1][1] / above_threshold_read_list[pos2][1]
                    mismatch_threshold = 2 #requires at least two mismatches 
                    if read_count_ratio  > 50:
                        mismatch_threshold = 3 #larger imbalances between reads changes the threshold for filtering out reads
                    if read_count_ratio  > 2500:
                        mismatch_threshold = 4
                    mismatches = count_mismatches(sorted_randomized_list[pos1][0], sorted_randomized_list[pos2][0], mismatch_threshold)
                    if mismatches < mismatch_threshold:
                        unique_dict[pos2] = 0
    for pos in range(len(above_threshold_read_list)):
        if unique_dict[pos] == 1:
            filtered_read_list.append(above_threshold_read_list[pos])
        

    print('%d filtered reads' %len(filtered_read_list))
    #print('%5.4f seconds for filtering' %(filter_end - filter_start))
    peptide_dict = {} #list of each sequence read encoding a given peptide sequence- indexed by peptide sequence
    ORF = 0
    for item in filtered_read_list:
        peptide = translate(item[0])
        if peptide.find('.') == -1:
            ORF += 1
            if peptide not in peptide_dict:
                peptide_dict[peptide] = [(item[0], item[1])]
            else:
                peptide_dict[peptide].append((item[0], item[1]))
    peptide_list = list(peptide_dict.keys())
    print('%d unique peptides' %len(peptide_list))
    peptide_list_passing_readcount_threshold = []
    for peptide in peptide_list:
        if peptide_dict[peptide][0][1] >= read_count_threshold:
            peptide_list_passing_readcount_threshold.append(peptide)
    pct_ORF = 100.0 * float(ORF) / float(len(filtered_read_list))
    print('%s\t%d\t%d\t%d\t%d\t%d\t%5.2f' %(input_filename, total, constant_region_correct, len(sorted_randomized_list), len(filtered_read_list), ORF, pct_ORF))
    input_file.close()
    output_list = []
    output_peptide_dict = {}
    for item in filtered_read_list:
        peptide = translate(item[0])
        if peptide.find('.') == -1:
            if peptide not in output_peptide_dict: #using dictionary to quickly filter out duplicates- not actually looking up info
                output_list.append((item[0], item[1], peptide))
                output_peptide_dict[peptide] = 1
                                   
    sorted_output_list = sorted(output_list, key=lambda x:x[1], reverse=True)
    sorted_peptide_list = []
    for item in sorted_output_list:
        output_line = '%s\t%d\t%s\n' %(item[0], item[1], item[2])
        output_all_file.write(output_line)
        sorted_peptide_list.append(item[2])
    output_all_file.close()

    truncated_peptide_list =  sorted_peptide_list[:number_of_peptides_to_plot]       
    
