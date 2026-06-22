
import os
import subprocess
import argparse
import shutil

############################################
"""
This script filters for the start of the trimmed R2 reads that at this point in the pipeline should correspond to the 
dinucleotide position. Filters for the correct donor dinucleotide used in the experiment, which is supplied as an argument. 
Found this step has helped remove some read artifacts:
1. Define the input directory as 'trimmed_fastq', move to 'trimmed_fastq_prior', 
    and create a new 'trimmed_fastq' directory as the output directory.
2. Perform cutadapt filtering on R2 files that should have the dinucleotide at the start of the read after previous 
    trimming steps. Perform separately for plus and minus files, minus files will have the reverse complement of the
    intended dinucleotide.
3. Move over the R1 files to the new trimmed_fastq directory.
4. Print a message to indicate that filtering is complete.
"""
############################################

#Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Path to input directory")
parser.add_argument("dinucleotide", help='Donor dinucleotide used in experiment')
args = parser.parse_args()

#Define directories
trimmed_fastq_prior = "trimmed_fastq_prior"
trimmed_fastq_dir = "trimmed_fastq"

#Rename existing trimmed_fastq to trimmed_fastq_prior and make new trimmed_fastq folder to seamlessy work with rest of pipeline
try:
    os.rename(trimmed_fastq_dir, trimmed_fastq_prior)
except FileNotFoundError:
    print(f"Error: '{trimmed_fastq_dir}' directory not found")
    exit(1)
except PermissionError:
    print(f"Error: insufficient permissions to rename '{trimmed_fastq_dir}'")
    exit(1)

os.makedirs(trimmed_fastq_dir, exist_ok=True)

#Get reverse complement of dinucleotide for minus files
complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
rev_comp = ''.join(complement[base] for base in reversed(args.dinucleotide.upper()))
dinucleotide = args.dinucleotide.upper()


#Filter R2 files for correct dinucleotide at start of read, retained reads will have lower-case dinucleotide at 
#beginning of reads
for file in os.listdir(trimmed_fastq_prior):
    if not (file.endswith(".fastq.gz") and "_R2_" in file):
        continue

    input_file = os.path.join(trimmed_fastq_prior, file)
    output_file = os.path.join(trimmed_fastq_dir, file)

    if "plus" in file:
        anchor_seq = dinucleotide
    elif "minus" in file:
        anchor_seq = rev_comp
    else:
        print(f"Skipping {file}: filename contains neither 'plus' nor 'minus'")
        continue

    try:
        subprocess.run([
            "cutadapt", input_file,
            "--cores", "16",
            "-g", f"^{anchor_seq}",
            "--action=lowercase",
            "--discard-untrimmed",
            "-o", output_file
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {file}: {e}")

#Copy R1 files to new trimmed_fastq directory
for file in os.listdir(trimmed_fastq_prior):
    if file.endswith(".fastq.gz") and "_R1_" in file:
        shutil.copy2(os.path.join(trimmed_fastq_prior, file), os.path.join(trimmed_fastq_dir, file))

print("Dinucleotide filtering complete")