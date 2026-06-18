from Bio import SeqIO
from multiprocessing import Pool
import pandas as pd
import argparse
import os
import gzip
from Levenshtein import distance

####################################################################
""" The code below does the following:
1. Creates a Pool object with 16 processes
2. Creates a list of tuples, containing the input and output files
3. Calls the filter_fastq function on each file, in parallel
    * keep all reads which contains sequence1 (start of read) but not sequence2 (unintegrated plasmid)
    * performs this separately for plus and minus files 
        (looks separately for corresponding top and bottom strands of AttP sequences)
    * allows 2 mismatches to the 5' sequence (23 nt) in search
    * allows 1 mismatch to the 3' sequence (15 nt) in search
4. Saves the reads that do not have both sequence1 and sequence2 to the output file
5. Saves the results of filter_fastq to a dataframe 
    (this is to track total number of reads and number of reads left after the plasmid filtering step)
6. Saves the dataframe to a csv file 
7. Print a message to indicate filtering is complete"""
####################################################################


#Define the two sequences to search for with sequence 1 being known and part of R2 and sequence 2 being the first sequence 
#after the dinucleotide corresponding to unintegrated plasmid. Sequences should have sequence 1 but throw out those with 
#sequence 2. Minus is the minus reaction version of the sequences (reads starting from the integrated plasmid, 
#downstream of the dinucleotide).
sequence1 = "GGTTTGTCTGGTCAACCACCGCG"
sequence2 = "CTCAGTGGTGTACGG"
sequence_minus_1 = "GGTTTGTACCGTACACCACTGAG"
sequence_minus_2 = "CGCGGTGGTTGACCA"

#match a string to expected pattern, allowing a max number of mismatches
def has_approximate_match(pattern, sequence, max_mismatches):
    window_size = len(pattern)
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        if distance(pattern, window) <= max_mismatches:
            return True
    return False

#Keep reads with seq1 but not seq2, allowing appropriate number of mismatches for pattern searched for
def filter_fastq(input_file, output_file, seq1, seq2):
    total_reads = 0
    filtered_reads = 0

    with gzip.open(input_file, "rt") as handle, gzip.open(output_file, "wt") as output_handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq_record = str(record.seq)
            if has_approximate_match(seq1, seq_record, 2) and not has_approximate_match(seq2, seq_record, 1):
                SeqIO.write(record, output_handle, "fastq")
                filtered_reads += 1

    return total_reads, filtered_reads

#Process R2 fastq.gz files from input directory, process separately for plus and minus files due to different string patterns
#expected. Look for "plus" and "minus" in file names for processing. Count reads making it through filter, this gives you the
#rough amount of starting usable reads that don't consist of unintegrated plasmid. Write out to attpFULL_filt_dir folder.
def main():
    #Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Path to input directory")
    args = parser.parse_args()

    attpFULL_filt_dir = os.path.join(os.getcwd(), "attpFULL_filt_dir")
    os.makedirs(attpFULL_filt_dir, exist_ok=True)

    num_processes = 16
    pool = Pool(processes=num_processes)
    counts_dict = {}
    
    results = []
    for files in os.listdir(args.input_dir):
        if "plus" in files and files.endswith("R2.fastq.gz"):
            input_file = os.path.join(args.input_dir, files)
            basename = os.path.splitext(os.path.basename(files))[0].replace("_R2.fastq", "")
            output_file = os.path.join(attpFULL_filt_dir, basename + "_attpFULLfilt_R2.fastq.gz")
            print(f"Processing file {input_file}")
            result = pool.apply_async(filter_fastq, (input_file, output_file, sequence1, sequence2))
            results.append((files, result))
            
        elif "minus" in files and files.endswith("R2.fastq.gz"):
            input_file = os.path.join(args.input_dir, files)
            basename = os.path.splitext(os.path.basename(files))[0].replace("_R2.fastq", "")
            output_file = os.path.join(attpFULL_filt_dir, basename + "_attpFULLfilt_R2.fastq.gz")
            print(f"Processing file {input_file}")
            result = pool.apply_async(filter_fastq, (input_file, output_file, sequence_minus_1, sequence_minus_2))
            results.append((files, result))

    pool.close()
    pool.join()

    
    for file, result in results:
        try:
            total, filtered = result.get()
            counts_dict[file] = {"total_reads": total, "filtered_reads": filtered}
        except Exception as e:
            print(f"Error processing {file}: {e}")


    df = pd.DataFrame.from_dict(counts_dict, orient='index')
    df.index = df.index.str.replace("_R2.fastq.gz", "", regex=False)
    df = df.reset_index()
    df.rename(columns={'index': 'file_name'}, inplace=True)
    df.to_csv("read_counts.csv", index=False)

    print("Plasmid filtering complete")



if __name__ == '__main__':
    main()
