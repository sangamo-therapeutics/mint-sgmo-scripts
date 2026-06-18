import os
import pybedtools
from pybedtools import BedTool

##############################################################################
"""
This script processes BAM files containing aligned sequencing reads and converts them to BED format

The code does the folowing:
Set input and output directory paths as "aligned" and "bed_files", respectively.
Iterate through files in the input directory & do the following:
    a. Check if the file ends with "_sorted.bam".
    b. Remove "_sorted" from the file's base name and define the output file path.
    c. Create a BedTool object from the input BAM file using pybedtools.BedTool().
    d. Convert the BAM file to a BED file using the bam_to_bed() method.
    e. Sort the BED file using sort(), and save it with saveas().
    f. Print a message to indicate bed files have been created and summarized.
"""
##############################################################################

#Define directories
input_dir = "aligned"
output_dir = "bed_files"

os.makedirs(output_dir, exist_ok=True)

#Convert bam files sorted by coordinates to bed files
for file in os.listdir(input_dir):
    if file.endswith("_sorted_name.bam"):
        print(f"processing: {file}")
        base = os.path.splitext(file)[0].replace("_sorted_name", "")
        try:
            bam_file = pybedtools.BedTool(os.path.join(input_dir, file))
            bed_file = bam_file.bam_to_bed(bedpe=True, mate1=True)
            bed_file.sort().saveas(os.path.join(output_dir, f"{base}.bed"))
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue

print("Bed files created")