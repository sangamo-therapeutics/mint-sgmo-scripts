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

#### usage: python Figure_1F_process_and_plot.py oligo_NNN_library_3-6-2024.txt Figure_1F_SATALKR_fq_files.txt

## Extended Data Figure 2B:
To reproduce the Bxb1 helix selection results in Extended Data Figure 2B you will need the follwing data files:
* ExtD_FIGURE_2B_GAC_batch2.fq (GAC)
* ExtD_FIGURE_2B_GTC_batch1.fq (GTC)
* ExtD_FIGURE_2B_CAG_batch4b.fq (CAG) [not uploaded to SRA yet]
* ExtD_FIGURE_2B_AGG_batch2.fq (AGG)
* ExtD_FIGURE_2B_AAT_batch4.fq (AAT)
* ExtD_FIGURE_2B_CTC_batch2.fq (CTC)
* ExtD_FIGURE_2B_CAA_batch1.fq (CAA)
* ExtD_FIGURE_2B_CAT_batch3.fq (CAT)
* ExtD_FIGURE_2B_CGC_batch4.fq (CGC)
* ExtD_FIGURE_2B_ACA_batch2.fq (ACA)

Generating the plots is a two step process where:
* the first step converts the raw sequence data into a sorted list of selected helix sequences
* the second step identifies enriched 4 residue patterns and then generates a plot from the 100 selected helix sequences that best match the enriched patterns

#### usage: python helix_selection_processing.py helix_library_template.txt ExtD_FIGURE_2B_AAT_batch4.fq
#### usage: python helix_pattern_finder_and_plot.py ExtD_FIGURE_2B_AAT_batch4_peptides_all.txt
  

## Extended Data Figure 2D:
To reproduce the Bxb1 hairpin selection results in Extended Data Figure 2E you will need the follwing data files:
* ExtD_FIGURE_2D_chr1_25477444L_GCCCCTTC_batch4.fq (GCCCCTTC)
* ExtD_FIGURE_2D_chr19_48971022L_GGGATTCC_batch3.fq (GGGATTCC)
* ExtD_FIGURE_2D_AAVS15032L_CTGAGCGC_batch3.fq (CTGAGCGC)
* ExtD_FIGURE_2D_AAVS15032R_GGGTTTGA_batch3.fq (GGGTTTGA)
* ExtD_FIGURE_2D_chr19_48971022R_TGCCTTCC_batch4.fq (TGCCTTCC)

Generating the plots is a multistep process where:
* the first step converts the raw sequence data into a sorted list of selected peptide sequences
* the second step identifies enriched 4 residue motifs for each of G, A, R, and P fixed at position 322
* the third step clusters the enriched patterns into separate motifs and generates a plot from the selected sequences that best match each motif

Note that some files contain multiple motifs and a separate plot is generated for each motif. 
In Extended Data Figure 2E the plot for GCCCCTTC represents motif3, the plot for TGCCTTCC represents motif4, and the other panels represent motif1 from their respective sample

#### usage: python hairpin_selection_processing.py hairpin_library_template.txt ExtD_FIGURE_2D_AAVS15032L_CTGAGCGC_batch3.fq
#### usage: python hairpin_pattern_finder.py ExtD_FIGURE_2D_AAVS15032L_CTGAGCGC_batch3_peptides_all.txt
#### usage: python hairpin_pattern_sort_and_plot.py ExtD_FIGURE_2D_AAVS15032L_CTGAGCGC_batch3_peptides_all_4res_patterns.txt

## Extended Data Figure 2F:
To reproduce the DNA specificity motif for the YRGSLP loop in Extended Data Figure 2F, you will need the following data files:
* ExtD_FIGURE_2F_YRGSLP_rep1.fq
* ExtD_FIGURE_2F_YRGSLP_rep2.fq
* ExtD_FIGURE_2F_YRGSLP_rep3.fq
* ExtD_FIGURE_2F_YRGSLP_rep4.fq

#### usage: python ExtD_Figure_2F_process_and_plot.py oligo_NN_library_loop.txt ExtD_FIGURE_2F_YRGSLP_fq_files.txt


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

#### usage: python ExtD_Figure_2G_process_and_plot.py oligo_NNNNNN_hairpin_library.txt ExtD_FIGURE_2G_FAGGGRKHPRYR_fq_files.txt
