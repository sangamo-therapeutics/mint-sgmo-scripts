#!/bin/bash

#This iteration of the code is used for processing potential genome-wide translocation sites with the TRAC on-target
#and with two replicates. To be used for processing data associated with SI tables 30-36.

set -eo pipefail

#VARIABLES

data_dir=$1

#Will look in directory for fastq.gz files and replaces suffix with new intermediate file format names, 
#saves R1 and R2 to separate folders for ease of processing.

#Execute this script using command: bash pipeline_genomic_TRAC_2_reps.sh .

#trim adapters and save R1 and R2 in separate folders
mkdir -p $data_dir/trimmed_fastq_R2

for file in $data_dir/*R2.fastq.gz; do
  filename=$(basename "$file" .fastq.gz)
  report_file="$data_dir/trimmed_fastq_R2/${filename}_trimming_report.txt"
  cutadapt -a AGATCGGAAGAGCGT --cores 16 -o "$data_dir/trimmed_fastq_R2/${filename}_trimmed.fastq.gz" "$file" --report minimal > "$report_file"
done

mkdir -p $data_dir/trimmed_fastq_R1

for file in $data_dir/*R1.fastq.gz; do
  filename=$(basename "$file" .fastq.gz)
  report_file="$data_dir/trimmed_fastq_R1/${filename}_trimming_report.txt"
  cutadapt -a AGATCGGAAGAGCGT --cores 16 -o "$data_dir/trimmed_fastq_R1/${filename}_trimmed.fastq.gz" "$file" --report minimal > "$report_file"
done

#Filter out expected sequence before dinucleotide in R2, output to genome_filt
python filter_genomic_reads_TRAC.py

mkdir -p $data_dir/genome_trimmed_fastq

#trim R2 files up to dinucleotide using expected genomic sequence, same sequence as previously filtered for in prior step

for file in $data_dir/genome_filt/*plus*.fastq.gz; do
  filename=$(basename "$file" .fastq.gz)
  report_file="$data_dir/genome_trimmed_fastq/${filename}_trimming_report.txt"
  cutadapt -g AATGGTGTCCAGGAGCCGAG --minimum-length 30 --cores 16 -o "$data_dir/genome_trimmed_fastq/${filename}.fastq.gz" "$file" --report minimal > "$report_file"
done

for file in $data_dir/genome_filt/*minus*.fastq.gz; do
  filename=$(basename "$file" .fastq.gz)
  report_file="$data_dir/genome_trimmed_fastq/${filename}_trimming_report.txt"
  cutadapt -g CTGGCCCTGGCAGGACCGAT --minimum-length 30 --cores 16 -o "$data_dir/genome_trimmed_fastq/${filename}.fastq.gz" "$file" --report minimal > "$report_file"
done

#remove plasmid reads in R2 using expected sequence after dinucleotide

python filter_plasmid_reads_genomic.py

python pair_align_reads_genomic.py

#follow the rest of the standard pipeline

python bed_formatting_paired.py

python R1dedupe_paired.py

python process_integration_loci_combine_strands.py

python plus_minus_2_rep_processing_one_file.py 

python update_coord.py


exit

