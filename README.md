# MINT scripts
This repository contains python scripts and configuration files to regenerate plots and analyses from Fauser et al. Nature Biotechnology 2026.

The raw data files can be found in the NCBI Sequence Read Archive under BioProject accession number PRJNA1450726

The motif plots use logomaker

## Figure 1F:
To reproduce the DNA specificity motif for the SATALKR helix in Figure 1F, you will need the following data files: 
* Figure_1F_SATALKR_rep1.fq
* Figure_1F_SATALKR_rep2.fq
* Figure_1F_SATALKR_rep3.fq
* Figure_1F_SATALKR_rep4.fq

###### usage: python Figure_1F_process_and_plot.py oligo_NNN_library_3-6-2024.txt Figure_1F_SATALKR_fq_files.txt


## Extended Data Figure 2F:
To reproduce the DNA specificity motif for the YRGSLP loop in Extended Data Figure 2F, you will need the following data files:
* ExtD_FIGURE_2F_YRGSLP_rep1.fq
* ExtD_FIGURE_2F_YRGSLP_rep2.fq
* ExtD_FIGURE_2F_YRGSLP_rep3.fq
* ExtD_FIGURE_2F_YRGSLP_rep4.fq

###### usage: python ExtD_Figure_2F_process_and_plot.py oligo_NN_library_loop.txt ExtD_FIGURE_2F_YRGSLP_fq_files.txt


## Extended Data Figure 2G:
To reproduce the DNA specificity motif for the FAGGG loop in Extended Data Figure 2F, you will need the following data files:
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep1_biorep1.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep1_biorep2.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep1_biorep3.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep1_biorep4.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep2_biorep1.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep2_biorep2.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep2_biorep3.fq
* ExtD_FIGURE_2G_FAGGGRKHPRYR_techrep2_biorep4.fq

###### usage: python ExtD_Figure_2G_process_and_plot.py oligo_NNNNNN_hairpin_library.txt ExtD_FIGURE_2G_FAGGGRKHPRYR_fq_files.txt
