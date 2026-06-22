import os
import pandas as pd
from pybedtools import BedTool

########################################################################
"""
This script combines both bed files of integration coordinates across both associated plus and minus reactions per sample, 
merges entries across both files that are within 50 bp of one another while maintaining the information on source file 
of integration events and alignment strands, then filters for at least 2 incidents of these events that were merged. 
This variation is written for mixed dinucleotide donors and no replicates.

The code does the following:
    a. Define the input and output directories as "bed_files_summed" and "bed_files_combined" respectively, 
    make the output directory if it doesn't exist
    b. In the input directory look for the 2 types of bed files associated with each TC sample ("plus" and "minus") 
    by searching for string descriptors in filenames, for files to be processed together these string descriptors can be
    the only difference in filenames
    c. Perform the combine_bedfiles function on each set of plus and minus files for a sample:
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

#Read in plus and minus files for 1 replicate, append file name to each file, combine all files, then merge entries within
#50 bp of one another, retain only locations with at least 2 instances in dataset
def combine_bedfiles(plus_filename, minus_filename):
    plus_file_path = os.path.join(input_directory, plus_filename)
    minus_file_path = os.path.join(input_directory, minus_filename)
    output_file_path = os.path.join(output_directory, plus_filename.replace("_summarized.bed", "_combined.bed"))

    #read bed files in as dataframes, add string descriptor of file source for each integration site
    plus_df = pd.read_csv(plus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    plus_df['source'] = 'plus'
    plus_bed = BedTool.from_dataframe(plus_df)
    minus_df = pd.read_csv(minus_file_path, sep='\t', names=['chr', 'start', 'stop', 'strand', 'count'])
    minus_df['source'] = 'minus'
    minus_bed = BedTool.from_dataframe(minus_df)

    #combine bed files, merge if within 50 bp, combine file sources in new column, filter for only sites with at least 
    #2 instances in dataset
    cat_bed = plus_bed.cat(*[minus_bed], postmerge=False)
    cat_bed_sort = cat_bed.sort()
    merged_bed = cat_bed_sort.merge(d=50, c=(6,4,2,2,5,5), o=('collapse','collapse','collapse','mean','collapse','mean'), delim = '_')
    filtered_bed = merged_bed.filter(lambda feature: '_' in feature.name)
    filtered_bed.saveas(output_file_path)


#Define directories
input_directory = 'bed_files_summed'
output_directory = 'bed_files_combined'

os.makedirs(output_directory, exist_ok=True)

#Get a list of plus and minus BED files
plus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("plus" in filename)]
minus_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.bed') and ("minus" in filename)]

#Check if plus and minus files associated with a TC sample are present and then process both to make a combined final bed file
for plus_filename in plus_files:
    minus_filename = plus_filename.replace("plus", "minus")
    if minus_filename in minus_files:
        combine_bedfiles(plus_filename, minus_filename)
    else:
        print(f"{plus_filename}: file is missing matching minus bed file")
    
print("Integration loci merged across plus/minus files, minimum 2 instances of integration event either on different strands or same strand in plus and minus files.")
