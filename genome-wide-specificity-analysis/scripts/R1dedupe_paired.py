import os
import pandas as pd

##############################################################################
"""
This script processes aligned paired read BED files, changes the read1 and read2 start to correspond to the start of 
the read or dinucleotide, respectively, irregardless of strand orientation, deduplicates reads based on unique R1 and 
R2 starts, then counts the number of R2 starts to make a summary BED file based on R2 information and number of unique 
R1 start locations (unique shear/ligation events).

The code does the folowing:
Set input and output directory paths as "bed_files" and "bed_files_R1dedupe", respectively.
Iterate through files in the input directory & do the following:
    a. Make output directory "bed_files_R1dedupe"
    b. In the input directory "bed_files" read in bed files
    c. Change R1 and R2 starts to correspond to start of read or dinucleotide, respectively, regardless of strand orientations
	d. Deduplicate on unique R1 and R2 starts to identify unique shear/ligation events in reactions
	e. Count the number of deduplicated instances per R2 starts (more quantitative representation of integration 
	strength than just read number)
	f. Make a new bed file with this summary and using R2 information (start, start + 1, strand)
	g. Write bed file to output directory with suffix "_R1dedupe.bed"
"""
##############################################################################

#Define directories
input_directory = "bed_files"
output_directory = "bed_files_R1dedupe"

os.makedirs(output_directory, exist_ok=True)

#Deduplicate reads according to R1 start location
def R1_dedupe(read_filename):
    read_file_path = os.path.join(input_directory, read_filename)

    #Extract the filename before ".bed" and make the output filename
    output_filename = read_filename.replace(".bed", "_R1dedupe.bed")
    output_file_path = os.path.join(output_directory, output_filename)

    #Read bed files in as dataframes
    read_df = pd.read_csv(read_file_path, sep='\t', names=['chr_R1', 'start_R1', 'stop_R1', 'chr_R2', 'start_R2', 'stop_R2', 'readid', 'qscore', 'strand_R1', 'strand_R2'])
    
    #Should be the case, but make sure chr_R1 == chr_R2
    read_df = read_df[read_df['chr_R1'] == read_df['chr_R2']]

    #Change the start location of R2 and R1 depending on strand location (negative uses end coordinates) 
    #to correspond to start of dinucleotide or start of read
    read_df['start_R1'] = read_df.apply(lambda x: x['stop_R1'] if x['strand_R1'] == '-' else x['start_R1'], axis=1)
    read_df['start_R2'] = read_df.apply(lambda x: x['stop_R2'] if x['strand_R2'] == '-' else x['start_R2'], axis=1)

    #Make new df with all required columns, deduplicate on R1/R2 starts, 
    #also need to keep R2 strand information for downstream processing
    dedupe_df = read_df[['chr_R1', 'start_R2', 'start_R1', 'strand_R2']].drop_duplicates()

    #Group by 'chr', and 'start_R2', and count to give the deduplicated numbers
    count_df = dedupe_df.groupby(['chr_R1', 'start_R2', 'strand_R2']).size().reset_index(name='count')

    #Create a new DataFrame with 'chr', 'start_R2', 'strand_R2', and 'count'
    new_df = count_df[['chr_R1', 'start_R2', 'strand_R2', 'count']].copy()

    #Increment the 'start_R2' values to get 'start_R2 + 1' as the new bed end coordinate that's just a dummy column to make 
    #bed format
    new_df['start_R2 + 1'] = new_df['start_R2'] + 1

    #Reorder the columns into bed format
    new_df = new_df[['chr_R1', 'start_R2', 'start_R2 + 1', 'strand_R2', 'count']]

    #Save the dataframe to output_file_path
    new_df.to_csv(output_file_path, sep='\t', header=False, index=False)


#Iterate through BED files in input directory
read_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed')]

for read_filename in read_files:
    R1_dedupe(read_filename)


print("R2 reads deduped on start of R1")


