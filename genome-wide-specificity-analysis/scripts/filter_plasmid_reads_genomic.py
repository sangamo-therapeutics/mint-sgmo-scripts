from Bio import SeqIO
from multiprocessing import Pool
import pandas as pd
import os
import gzip
from Levenshtein import distance

####################################################################
""" The code below does the following:
1. Creates a Pool object with 16 processes
2. Creates a list of tuples, containing the input and output files
3. Calls the filter_fastq function on each file, in parallel
    * removes all reads which contains either sequence1 or sequence2 (plasmid)
    * allows 1 mismatch to the 15 nt sequence in search
4. Saves the reads that do not have either sequence1 and sequence2 to the output file
5. Saves the results of filter_fastq to a dataframe 
    (this is to track total number of reads and number of reads left after the filtering step)
6. Saves the dataframe to a csv file 
7. Print a message to indicate filtering is complete"""
####################################################################


#Define the two sequences to search for (either the 5' end, top strand or 3' end, bottom strand of the plasmid once integrated)
plasmid_3_AttP = 'CTCAGTGGTGTACGG'
plasmid_5_AttP = 'CGCGGTGGTTGACCA'


#match a string to expected pattern, allowing a max number of mismatches
def has_approximate_match(pattern, sequence, max_mismatches):
    window_size = len(pattern)
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        if distance(pattern, window) <= max_mismatches:
            return True
    return False


#Keep reads without a match to seq1 or seq2, allowing appropriate number of mismatches for pattern searched for
def filter_fastq(input_file, output_file, seq1, seq2):
    total_reads = 0
    filtered_reads = 0

    with gzip.open(input_file, "rt") as handle, gzip.open(output_file, "wt") as output_handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq_record = str(record.seq)
            if not has_approximate_match(seq1, seq_record, 1) and not has_approximate_match(seq2, seq_record, 1):
                SeqIO.write(record, output_handle, "fastq")
                filtered_reads += 1

    return total_reads, filtered_reads


#Process R2 fastq.gz files from genome_trimmed_fastq. Count reads making it through filter, this gives you the rough amount of 
#starting usable reads that don't consist of unintegrated plasmid. Write out to plasmid_filt folder.
def main():
    plasmid_filt_dir = os.path.join(os.getcwd(), "plasmid_filt")
    os.makedirs(plasmid_filt_dir, exist_ok=True)
    genome_filt_dir = "genome_trimmed_fastq"

    num_processes = 16
    pool = Pool(processes=num_processes)
    counts_dict = {}
    
    results = []

    files_list = [f for f in os.listdir(genome_filt_dir) if f.endswith('.fastq.gz')]

    for files in files_list:
        input_file = os.path.join(genome_filt_dir, files)
        print(f"Processing file: {input_file}")
        basename = os.path.basename(files).replace(".fastq.gz", "")
        output_file = os.path.join(plasmid_filt_dir, basename + ".fastq.gz")
        result = pool.apply_async(filter_fastq, (input_file, output_file, plasmid_3_AttP, plasmid_5_AttP))
        results.append((files, result))

    pool.close()
    pool.join()

    
    for files, result in results:
        try:
            total, filtered = result.get()
            counts_dict[files] = {"total_reads": total, "filtered_reads": filtered}
            print(f"Processed file {files}, total reads: {total}, filtered reads: {filtered}")
        except Exception as e:
            print(f"Error processing {files}: {e}")
            

    df = pd.DataFrame.from_dict(counts_dict, orient='index')
    df.index = df.index.str.replace("_R2_trimmed.fastq.gz", "")
    df = df.reset_index()
    df.rename(columns={'index': 'file_name'}, inplace=True)
    df.to_csv(os.path.join(plasmid_filt_dir, 'read_counts.csv'), index=False)

    print("Filtering complete!")



if __name__ == '__main__':
    main()
