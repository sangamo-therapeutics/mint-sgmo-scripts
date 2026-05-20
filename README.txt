This repository contains python scripts and configuration files to regenerate plots and analyses from Fauser et al. Nature Biotechnology 2026.

The raw data files can be found in the NCBI Sequence Read Archive under BioProject accession number PRJNA1450726

Figure 1F:
To reproduce the DNA specificity motif for the SATALKR helix in Figure 1F, you will need the following data files: 
Figure_1F_SATALKR_rep1.fq
Figure_1F_SATALKR_rep2.fq
Figure_1F_SATALKR_rep3.fq
Figure_1F_SATALKR_rep4.fq

The script can be run with the following command:
python Figure_1F_process_and_plot.py oligo_NNN_library_3-6-2024.txt Figure_1F_SATALKR_fq_files.txt


