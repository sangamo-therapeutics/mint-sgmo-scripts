import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor


########################################################################
"""
This script processes the previously made bed files with deduplicated and summed number of reads, 
separated entries per strand (necessary for mixed dinucleotide reactions where one PCR reaction can
amplify both a top and bottom strand integration event based on the top or bottom strand dinucleotide),
groups entries that are within 50 bp of one another, chooses the read start position 
with the most reads in a group, and associates that entry with summed reads per group, 
combines the results for top and bottom strand alignments and saves the results in a new sorted BED file

The code does the following:
    a. Set input and output directory paths as 'bed_files_R1dedupe' and 'bed_files_summed', respectively.
    b. Create the output directory if it doesn't exist using os.makedirs().
    c. In the input directory, read in "_R1dedupe.bed" files
    d. Perform the process_bedfile() function to process each BED file:
        a. Read the BED file into a pandas DataFrame.
        b. Separate BED file into two DataFrames according to alignment strand.
        c. Perform the process_strand_group() function to process each BED file:
            a. Sort DataFrames by chromosome, start position, and descending value.
            b. Group loci in a 50 bp window, choose the start position with the highest value, 
            sum values in a group to report with start position.
        d. Combine the top and bottom strand processed bed files, sort
        e. Save the resulting DataFrame to a new BED file with the suffix "_summarized.bed".
    e. Print a message to indicate bed files have been processed
"""
########################################################################

#Group loci within a 50 bp window, choose the start position with the highest value to report but sum all values in a group. 
def process_strand_group(df):
    #Sort by chromosome, start position, and descending value
    df_sorted = df.sort_values(by=['chr', 'start', 'value'], ascending=[True, True, False]).copy()

    #Group loci within a 50 bp window, choose the start position with the highest value to report but sum all values in a group.
    #This worked well with raw read numbers but with R1 dedupe the values are lower and less varied so there are instances 
    #where there is no location within a 50 bp window with the highest value to report. In this scenario, 
    #the last chromosomal location of the group in the list is reported.
    result_rows = []
    curr_window = []
    for _, row in df_sorted.iterrows():
        if not curr_window or row['chr'] == curr_window[-1]['chr'] and row['start'] - curr_window[-1]['start'] <= 50:
            curr_window.append(row)
        else:
            window_sum = sum(r['value'] for r in curr_window)
            max_row = max(curr_window, key=lambda x: x['value'])
            max_row['value'] = window_sum
            result_rows.append(max_row)
            curr_window = [row]

    #Add the last window
    if curr_window:
        window_sum = sum(r['value'] for r in curr_window)
        max_row = max(curr_window, key=lambda x: x['value'])
        max_row['value'] = window_sum
        result_rows.append(max_row)

    return pd.DataFrame(result_rows)

#Process bed files separately according to R2 strand
def process_bedfile(filename):
    #Read in files in input directory, define corresponding output files
    file_path = os.path.join(input_directory, filename)
    output_file_path = os.path.join(output_directory, filename.replace("_R1dedupe.bed", "_summarized.bed"))

    #Read the BED file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chr', 'start', 'end', 'strand', 'value'])

    #Split the DataFrame into two groups based on strand
    df_plus = df[df['strand'] == '+']
    df_minus = df[df['strand'] == '-']

    #Summarize reads in 50 bp window per strand
    df_plus_processed = process_strand_group(df_plus)
    df_minus_processed = process_strand_group(df_minus)

    #Combine the summarized DataFrames
    df_combined = pd.concat([df_plus_processed, df_minus_processed])

    #Sort the combined DataFrame by chromosome and start position in ascending order
    df_combined.sort_values(by=['chr', 'start'], ascending=[True, True], inplace=True)

    #Save the combined DataFrame to a new BED file
    df_combined.to_csv(output_file_path, sep='\t', header=False, index=False)

#Define input and output directories, make output directory
input_directory = "bed_files_R1dedupe"
output_directory = 'bed_files_summed'

os.makedirs(output_directory, exist_ok = True)

#Get a list of BED files in input directory
bed_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed')]

#Process BED files using multiple processors
with ProcessPoolExecutor() as executor:
    executor.map(process_bedfile, bed_files)

print("Integration loci defined, summed counts across 50 bp window per R2 strand")
