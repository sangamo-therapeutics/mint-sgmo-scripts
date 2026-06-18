#!/bin/bash

#This iteration of the code is used for processing potential genome-wide integration sites using a single dinucleotide donor 
#and with two replicates. To be used for processing data associated with SI tables 10-13 (CA dinucleotide donor) and 23-29
#(GT dinucleotide donor).

set -eo pipefail

#VARIABLES, directory where files are located and the donor dinucleotide used in the experiment (in this case either CA or GT)

data_dir=$1
dinucleotide=$2

#Will look in directory for *R2.fastq.gz files
#Execute this script using command: bash pipeline_nonpalindromic_nonmixed_2_reps.sh . <donor dinucleotide>

python filter_plasmid_reads.py $data_dir

python trimm_adapters_R1_R2.py $data_dir

python filter_nt.py $data_dir $dinucleotide

python pair_align_reads.py

python bed_formatting_paired.py

python R1dedupe_paired.py

python process_integration_loci_combine_strands.py

python plus_minus_2_rep_processing_one_file.py

python update_coord.py


exit
