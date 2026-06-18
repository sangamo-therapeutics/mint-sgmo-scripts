
import os
import subprocess
import argparse

############################################
"""
This script trims the start of the reads up to the dinucleotide to enable genome alignment:
1. Define the input directory as 'attpFULL_filt_dir', the output directory as 'trimmed_fastq', 
    and create the output directory if it doesn't exist.
2. Set the quality control (QC) report file path as 'trimmingQC.txt'.
3. Iterate through R2.fastq.gz files in the input directory:
   a. Print the processing file name.
   b. Extract the file base name and set the temporary and final output file paths.
   c. Removes up to the dinucleotide for plus and minus files via trimming the 5' AttP separately for each
   d. Trim adaptor sequences from the 3' end of each read and ensure a minimum read length of 30 after 5' and 3' 
    file (different sequences) and trimming using cutadapt, 
    while appending the QC report to the 'trimmingQC.txt' file.
4. Iterate through R1.fastq.gz files in the parent directory:
    a. Print the processing file name.
    b. Extract the file base name and set the temporary and final output file paths.
    c. Removes the beginning of R2 corresponding to the plasmid if present at the end of R1 (short reads) 
    for plus and minus files, otherwise this throws off alignment
    d. Trim adaptor sequences from the 3' end of each read, requiring 30 bp leftover for R1 alignment, 
    while appending the QC report to the 'trimmingQC.txt' file.
5. Remove temporary files from the output directory.
6. Print a message to indicate that trimming is complete.
"""
############################################

#Define directories and files
attpFULL_filt_dir = "attpFULL_filt_dir"
trimmed_fastq_dir = "trimmed_fastq"
os.makedirs(trimmed_fastq_dir, exist_ok=True)
qc_file = "trimmingQC.txt"


#Function to remove everything up to 5' AttP in R2
def R2_trimming(input_file, output_file, specific_seq, qc_file):
    #Remove everything up to the specific sequence
    subprocess.run([
        "cutadapt",
        "-g", specific_seq,
        "-o", output_file,
        "--cores", "16",
        input_file
    ], stdout=qc_file)


#Function to perform adapter trimming
def adapter_trimming(input_file, output_file, qc_file):
    #Trim adapter sequences from 3' end of each read
    subprocess.run([
        "cutadapt",
        "-a", "AGATCGGAAGAGCGT",  #Adapter sequence
        "--minimum-length", "30",  #Minimum length after trimming
        "-o", output_file,
        "--cores", "16",
        input_file,
        "--report", "minimal"
    ], stdout=qc_file)


#Function to remove plasmid at 3' of R1 if present (short reads, will throw off alignment)
def R1_trimming(input_file, output_file, specific_seq, qc_file):
    #Trim adapter sequences from 3' end of each read
    subprocess.run([
        "cutadapt",
        "-a", specific_seq,  #Beginning of plasmid seq after dinucleotide
        "-o", output_file,
        "--cores", "16",
        input_file,
        "--report", "minimal"
    ], stdout=qc_file)


#Function to process R2 files
def process_R2files(file, input_dir, trimmed_fastq_dir, specific_seq):
    base = file.replace("attpFULLfilt_R2.fastq.gz", "")
    temp_file = os.path.join(trimmed_fastq_dir, f"{base}_temp.fastq.gz")
    final_file = os.path.join(trimmed_fastq_dir, f"{base}R2_trimmed.fastq.gz")

    #Remove 5' up to dinucleotide in R2
    with open(qc_file, "a") as qc:
        R2_trimming(os.path.join(input_dir, file), temp_file, specific_seq, qc)

    #Perform 3' adapter trimming step
    with open(qc_file, "a") as qc:
        adapter_trimming(temp_file, final_file, qc)
    
    #Remove temporary file
    os.remove(temp_file)


#Function to process R1 files
def process_R1files(file, input_dir, trimmed_fastq_dir, specific_seq):
    base = file.replace(".fastq.gz", "")
    temp_file = os.path.join(trimmed_fastq_dir, f"{base}_temp.fastq.gz")
    final_file = os.path.join(trimmed_fastq_dir, f"{base}_trimmed.fastq.gz")

    #Perform 3' plasmid trimming step in R1
    with open(qc_file, "a") as qc:
        R1_trimming(os.path.join(input_dir, file), temp_file, specific_seq, qc)

    #Perform adapter trimming step
    with open(qc_file, "a") as qc:
        adapter_trimming(temp_file, final_file, qc)
    
    #Remove temporary file
    os.remove(temp_file)


#Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Path to input directory")
args = parser.parse_args()

#Iterate through the files in the attpFULL_filt_dir for R2 files
for file in os.listdir(attpFULL_filt_dir):
    if file.endswith("_R2.fastq.gz"):
        if "plus" in file:
            specific_seq = "GGTTTGTCTGGTCAACCACCGCG"
        elif "minus" in file:
            specific_seq = "GGTTTGTACCGTACACCACTGAG"
        else:
            print(f"Skipping {file}: filename contains neither 'plus' nor 'minus'")
            continue
        
        print(f"Processing {file}")
        process_R2files(file, attpFULL_filt_dir, trimmed_fastq_dir, specific_seq)

#Iterate through the files in the input directory for R1 files
for file in os.listdir(args.input_dir):
    if file.endswith("_R1.fastq.gz"):
        if "plus" in file:
            specific_seq = "CGCGGTGGTTGACCA"
        elif "minus" in file:
            specific_seq = "CTCAGTGGTGTACGG"
        else:
            print(f"Skipping {file}: filename contains neither 'plus' nor 'minus'")
            continue
        
        print(f"Processing {file}")
        process_R1files(file, args.input_dir, trimmed_fastq_dir, specific_seq)

print("Trimming complete")