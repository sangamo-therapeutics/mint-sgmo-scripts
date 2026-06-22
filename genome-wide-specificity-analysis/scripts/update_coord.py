import os
import glob
from pybedtools import BedTool

########################################################################
"""
This script reads in the output bed files from the specificity analysis and formats them to give 38 bp long putative attB 
off-target sequences. It does so by changing the start and end coordinates so that getfasta produces the right sequence with 
the dinucleotide centered when possible. For files that have reads in only 1 strand, if it's aligned to the bottom strand, 
the start coordinate is off by 1 and thus those are processed differently. For edge cases with atypical start and end sites 
compared to the expected pattern that may correspond to false positives, bad read quality/alignment/clipping, 
or have strange alignments based on the K562 genome, the given sequence will be difficult to make accurate without personal 
inspection into each.
The code does the following:
    a. Define the input and output directories as "bed_files_combined" and "bed_files_coord" respectively, 
    make the output directory if it doesn't exist
    b. In the input directory look for bed files
    c. Perform the rewrite coord function on each bed file:
        a. Read in each bed file as a bedtool bed file
        b. Rewrite start and end coordinates to center the theoretical dinucleotide in sequence, performed separately for 
        files that either have any + alignment or none due to the way bedtools outputs reads with - strand alignment 
        (reads that only align to the bottom strand that are off by one start coordinate)
        c. Getfasta for indicated coordinates
        e. Write out as a new bed file with the suffix "_coord" attached to the file name
    d. Prints a message to indicate all processing has been performed on samples.
"""
########################################################################

# Define directories
input_directory = "bed_files_combined"
output_directory = "bed_files_coord"
fasta_ref = "/seqdb/fasta/hg38/hg38all.fa"

# Column indices based on the previous output
IDX_CHR = 0
IDX_MERGED_START = 1
IDX_MERGED_END = 2
IDX_MERGED_SAMPLES = 3
IDX_MERGED_STRANDS = 4
IDX_STARTS = 5
IDX_MEAN_START = 6
IDX_COUNTS = 7
IDX_AVERAGE_COUNTS = 8

def rewrite_coordinates(in_bed, out_bed):
    """
    Rewrite coordinates based on strand and append fasta sequence as a new column.
    """
    # Load BED file with pybedtools
    bed = BedTool(in_bed)

    new_intervals = []
    for interval in bed:
        strand_info = interval[IDX_MERGED_STRANDS]

        start = int(interval[IDX_MERGED_START])
        # Current start coordinates correspond to putative dinucleotide positions
        # Compute new start based on strand, attB is 38 bp long
        # lack of any '+' is processed differently
        if "+" in strand_info:
            new_start = start + 1 - 19
        else:
            new_start = start - 1 - 19

        new_end = new_start + 38

        # Update interval fields
        fields = list(interval)
        fields[IDX_MERGED_START] = str(new_start)
        fields[IDX_MERGED_END] = str(new_end)

        # Use BedTool to get sequence from reference
        temp_bed = BedTool([fields])
        fasta_file = temp_bed.sequence(fi=fasta_ref, s=True)  # returns a BedTool object pointing to a FASTA file

        # Read sequence from fasta file
        with open(fasta_file.seqfn) as f:
            lines = f.readlines()
            # FASTA format: >header\nSEQUENCE\n
            sequence = "".join(line.strip() for line in lines if not line.startswith(">")).upper()

        # Append sequence as new column
        fields.append(sequence)
        new_intervals.append(fields)

    # Save all intervals with sequences
    BedTool(new_intervals).saveas(out_bed)


def main():
    os.makedirs(output_directory, exist_ok=True)

    bed_files = glob.glob(os.path.join(input_directory, "*.bed"))

    if not bed_files:
        raise RuntimeError(f"No .bed files found in {input_directory}")

    for bed_file in bed_files:
        base = os.path.basename(bed_file).rsplit(".bed", 1)[0]
        out_bed = os.path.join(output_directory, f"{base}_coord.bed")
    try:
        rewrite_coordinates(bed_file, out_bed)
        print(f"Processed: {bed_file} → {out_bed}")
    except Exception as e:
        print(f"Error processing {bed_file}: {e}")


if __name__ == "__main__":
    main()