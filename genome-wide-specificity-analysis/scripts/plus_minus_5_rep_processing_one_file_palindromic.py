import os
import pandas as pd
from pybedtools import BedTool

########################################################################
"""
This script combines all bed files of integration coordinates across all associated plus and minus reactions and replicates, 
merges entries across all files that are within 50 bp of one another while maintaining the information on source file 
of integration events and alignment strands, then filters for at least 2 incidents of these events that were merged. 
This variation is written for mixed dinucleotide donors and 5 replicates. 

The code does the following:
    a. Define the input and output directories as "bed_files_summed" and "bed_files_combined" respectively, 
    make the output directory if it doesn't exist
    b. In the input directory look for the 10 types of bed files associated with each sample conditions
    ("plus" and "minus" for each "_rep#_" file) by searching for string descriptors in filenames, 
    for files to be processed together these string descriptors can be the only difference in filenames
    c. Perform the combine_bedfiles function on each set of plus and minus and replicate files for each set of similar samples:
        a. Read in each bed file as a pandas df, add extra column to keep track of file source, then convert to bedtool bed file
        b. Combine and sort both bed file entries
        c. Merge all entries that are within 50 bp of one another, for merged entries concatenate file source, alignment strands,
        start locations, and counts with "_", for merged entries also take the average start location and reads
        d. Keep only entries merged from at least 2 files
        e. Write out as combined summary bed file with columns: "chr", "merged_start", "merged_end", 
        "files_merged", "strands_merged", "file_starts", "mean_start", "file_counts", "mean_counts"
    d. Prints a message to indicate all processing has been performed on samples.
"""
########################################################################

#Read in plus and minus files for 5 replicates, append file name to each file, combine all files, then merge entries within
#50 bp of one another, retain only locations with at least 2 instances in dataset
def combine_bedfiles(r1_plus_filename, r2_plus_filename, r1_minus_filename, r2_minus_filename, r3_plus_filename, 
    r3_minus_filename, r4_plus_filename, r4_minus_filename, r5_plus_filename, r5_minus_filename):
    r1_plus_file_path = os.path.join(input_directory, r1_plus_filename)
    r2_plus_file_path = os.path.join(input_directory, r2_plus_filename)
    r1_minus_file_path = os.path.join(input_directory, r1_minus_filename)
    r2_minus_file_path = os.path.join(input_directory, r2_minus_filename)
    r3_plus_file_path = os.path.join(input_directory, r3_plus_filename)
    r4_plus_file_path = os.path.join(input_directory, r4_plus_filename)
    r3_minus_file_path = os.path.join(input_directory, r3_minus_filename)
    r4_minus_file_path = os.path.join(input_directory, r4_minus_filename) 
    r5_plus_file_path = os.path.join(input_directory, r5_plus_filename)
    r5_minus_file_path = os.path.join(input_directory, r5_minus_filename)  
    output_file_path = os.path.join(output_directory, r1_plus_filename.replace("_summarized.bed", "_combined.bed"))

    #read 10 sets of bed files in as dataframes, add string descriptor of file source for each integration site
    r1_plus_df = pd.read_csv(r1_plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r1_plus_df['source'] = 'r1-plus'
    r1_plus_bed = BedTool.from_dataframe(r1_plus_df)
    r2_plus_df = pd.read_csv(r2_plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r2_plus_df['source'] = 'r2-plus'
    r2_plus_bed = BedTool.from_dataframe(r2_plus_df)
    r1_minus_df = pd.read_csv(r1_minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r1_minus_df['source'] = 'r1-minus'
    r1_minus_bed = BedTool.from_dataframe(r1_minus_df)
    r2_minus_df = pd.read_csv(r2_minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r2_minus_df['source'] = 'r2-minus'
    r2_minus_bed = BedTool.from_dataframe(r2_minus_df)
    r3_plus_df = pd.read_csv(r3_plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r3_plus_df['source'] = 'r3-plus'
    r3_plus_bed = BedTool.from_dataframe(r3_plus_df)
    r4_plus_df = pd.read_csv(r4_plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r4_plus_df['source'] = 'r4-plus'
    r4_plus_bed = BedTool.from_dataframe(r4_plus_df)
    r3_minus_df = pd.read_csv(r3_minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r3_minus_df['source'] = 'r3-minus'
    r3_minus_bed = BedTool.from_dataframe(r3_minus_df)
    r4_minus_df = pd.read_csv(r4_minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r4_minus_df['source'] = 'r4-minus'
    r4_minus_bed = BedTool.from_dataframe(r4_minus_df)
    r5_plus_df = pd.read_csv(r5_plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r5_plus_df['source'] = 'r5-plus'
    r5_plus_bed = BedTool.from_dataframe(r5_plus_df)
    r5_minus_df = pd.read_csv(r5_minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    r5_minus_df['source'] = 'r5-minus'
    r5_minus_bed = BedTool.from_dataframe(r5_minus_df)

    #combine 10 bed files, merge if within 50 bp, combine file sources in new column, filter for only sites with at least 
    #2 instances in dataset
    cat_bed = r1_plus_bed.cat(r2_plus_bed, r1_minus_bed, r2_minus_bed, r3_plus_bed, r3_minus_bed, r4_plus_bed, r4_minus_bed, 
        r5_plus_bed, r5_minus_bed, postmerge=False)
    cat_bed_sort = cat_bed.sort()
    merged_bed = cat_bed_sort.merge(d=50, c=(6,4,2,2,5,5), o=('collapse','collapse','collapse','mean','collapse','mean'), 
        delim = '_')
    filtered_bed = merged_bed.filter(lambda feature: '_' in feature.name)
    filtered_bed.saveas(output_file_path)

#Define directories
input_directory = 'bed_files_summed'
output_directory = 'bed_files_combined'

os.makedirs(output_directory, exist_ok=True)

#Get a list of plus and minus BED files for five replicates
r1_plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename) and ('_rep1_' in filename)]
r2_plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename) and ('_rep2_' in filename)]
r1_minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename) and ('_rep1_' in filename)]
r2_minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename) and ('_rep2_' in filename)]
r3_plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename) and ('_rep3_' in filename)]
r4_plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename) and ('_rep4_' in filename)]
r3_minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename) and ('_rep3_' in filename)]
r4_minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename) and ('_rep4_' in filename)]
r5_plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename) and ('_rep5_' in filename)]
r5_minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename) and ('_rep5_' in filename)]

#Check if plus and minus files across the five replicates associated with similar samples are present and then process all
#to make a combined final bed file
for r1_plus_filename in r1_plus_files:
    r2_plus_filename = r1_plus_filename.replace("_rep1_", "_rep2_")
    r3_plus_filename = r1_plus_filename.replace("_rep1_", "_rep3_")
    r4_plus_filename = r1_plus_filename.replace("_rep1_", "_rep4_")
    r5_plus_filename = r1_plus_filename.replace("_rep1_", "_rep5_")
    r1_minus_filename = r1_plus_filename.replace("plus", "minus")
    r2_minus_filename = r1_plus_filename.replace("_rep1_", "_rep2_").replace("plus", "minus")
    r3_minus_filename = r1_plus_filename.replace("_rep1_", "_rep3_").replace("plus", "minus")
    r4_minus_filename = r1_plus_filename.replace("_rep1_", "_rep4_").replace("plus", "minus")
    r5_minus_filename = r1_plus_filename.replace("_rep1_", "_rep5_").replace("plus", "minus")
    if r2_plus_filename in r2_plus_files and r3_plus_filename in r3_plus_files and r4_plus_filename in r4_plus_files and r5_plus_filename in r5_plus_files and r1_minus_filename in r1_minus_files and r2_minus_filename in r2_minus_files and r3_minus_filename in r3_minus_files and r4_minus_filename in r4_minus_files and r5_minus_filename in r5_minus_files:
        combine_bedfiles(r1_plus_filename, r2_plus_filename, r3_plus_filename, r4_plus_filename, r5_plus_filename, 
            r1_minus_filename, r2_minus_filename, r3_minus_filename, r4_minus_filename, r5_minus_filename)
    else:
        print(f"{r1_plus_filename}: file is missing a matching bed file in group")
    
print("Integration loci merged across plus/minus files for 5 replicates indicated by 'rep#' in names, minimum 2 instances of integration events.")