#!/bin/bash

#This iteration of the code is used for processing potential genome-wide integration sites using mixed dinucleotide donors 
#with only one replicate. To be used for processing data associated with SI tables 8 and 9.

set -eo pipefail

#VARIABLES

data_dir=$1

#Execute this script using command: bash pipeline_mixednt_1_rep.sh .
#Will look in directory for *R2.fastq.gz files

python filter_plasmid_reads.py $data_dir

python trimm_adapters_R1_R2.py $data_dir

python pair_align_reads.py

python bed_formatting_paired.py

python R1dedupe_paired.py

python process_integration_loci.py

python plus_minus_processing_one_file.py

python update_coord.py


exit
