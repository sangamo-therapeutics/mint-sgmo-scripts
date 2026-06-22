# README for processing all genome-wide specificity analysis code

The code here is required to generate all genome-wide specificity data in the paper (all data associated with Supplementary Tables 6, 8-13, 23-29, and 30-36).



## Necessary file properties for all analyses
There are 4 iterations of this code to account for experimental differences. All account for files with the following properties:
 - Each biological sample is analyzed using separate "**plus**" and "**minus**" primers in the assay that correspond to either direction of amplification with respect to the donor or on-target genomic integration site (depending on the assay). In these cases, "plus" refers to upstream of the dinucleotide and "minus" downstream for either the donor, for the genome-wide integration site nomination assay, or on-target genomic integration site, for the genome-wide translocation nomination assay.
 - If using replicates, each file has "**\_rep#\_**" in it's name to correspond to which one.
 - Files required for the analysis are provided in separate R1 and R2 format in the SRA and are required to be in unpaired format to run the analysis scripts. The files must end in "**R1.fastq.gz**" or "**R2.fastq.gz**".
 - Make sure all files associated with a set of experimental conditions (same reagents) are all named the same except for these variable names (plus/minus and rep#), separated by underscores, and do not have any other instance of these strings in the filenames. For example, here are all necessary files associated with Supplementary Table 10:
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**plus_rep1_R1.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**plus_rep1_R2.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**plus_rep2_R1.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**plus_rep2_R2.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**minus_rep1_R1.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**minus_rep1_R2.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**minus_rep2_R1.fastq.gz**
   - SI_TABLE_10_MACO1L_MACO1R_bxb1_**minus_rep2_R2.fastq.gz**

## Pipeline-data mapping
For genome-wide **integration** site nomination datasets that were generated using **mixed dinucleotide donors**, use `pipeline_mixednt_5_reps.sh` and `pipeline_mixednt_1_rep.sh`..
- For datasets with **5 replicates** use `pipeline_mixednt_5_reps.sh`. This is used on data associated with Supplementary Table 6.
- For datasets with only **1 replicate** use `pipeline_mixednt_1_reps.sh`. This is used on data associated with Supplementary Tables 8-9.

For genome-wide **integration** site nomination datasets that used a **single dinucleotide donor** and **2 replicates** use `pipeline_nonpalindromic_nonmixed_2_reps.sh`. Use it on data associated with Supplementary tables 10-13 and 23-29.

For genome-wide **translocation** site nomination datasets that used 2 replicates and correspond to interrogating potential recombinations at the TRAC target site use `pipeline_genomic_TRAC_2_reps.sh`. Use it on data associated with Supplementary Tables 30-36.

## Docker container set-up
The docker container set-up should be performed once and will allow the execution of all specificity scripts. Refer to the Docker README before executing.

## Processing files
1. Move read files associated with each experiment into their own directory as they will require different pipelines and called variables. Each script looks for R1 and R2.fastq.gz files in the directory it is run in and processes them all together.
   - **Experiment #1**: Supplementary Table 6 files
   - **Experiment #2**: Supplementary Tables 8-9 files
   - **Experiment #3**: Supplementary Tables 10-13 files
   - **Experiment #4**: Supplementary Tables 23-29 files
   - **Experiment #5**: Supplementary Tables 30-36 files
2. Move all pipelines and scripts into each experiment folder. If desired, you can also move only the pipeline and scripts necessary for each experiment but either way the pipeline executed will only use scripts necessary for running itself.
3. Execute the necessary pipeline command for each experiment:
   - Run `bash pipeline_mixednt_5_reps.sh .` for all files associated with Supplementary Table 6.
   - Run `bash pipeline_mixednt_1_reps.sh .` for all files associated with Supplementary Tables 8-9.
   - Run `bash pipeline_nonpalindromic_nonmixed_2_reps.sh . CA` for all files associated with Supplementary Tables 10-13.
   - Run `bash pipeline_nonpalindromic_nonmixed_2_reps.sh . GT` for all files associated with Supplementary Tables 23-29.
   - Run `bash pipeline_genomic_TRAC_2_reps.sh .` for all files associated with Supplementary Tables 30-36.
4. Final bed files are found in the bed_files_coord folder and should correspond to those in the Supplementary Tables, with additional columns of information. The final output is for the putative 38 bp attB site with the dinucleotide centered and supporting signal information from the assay. The columns are: **chromosome**, **start**, **end**, samples with signal (replicate # and plus/minus primer reaction), strands of the samples with signal, original start locations of the samples with signal (usually corresponds to before or after dinucleotide), average start location of the samples with signal, integration counts of the samples with signal, **average integration counts of the signal**, and the **putative attB site sequence** corresponding to the signal. All merged data columns have information separated by "_" and are combined in the same order per column. Bolded column descriptions correspond to the data retained in the published Supplementary Tables.
   - For example, this entry:
   > chr1	22229736	22229774	r3\-plus\_r4\-minus\_r1\-minus\_r3\-plus	+\_\-\_\-\_\-	22229754\_22229756\_22229756\_22229756	22229755.5	3\_1\_1\_1	1.5	CCTCTCTTTCTAACTGAACTGTCAGTTACACAAGCTCA
   - corresponds to a putative integration site with 1.5 average integration counts derived from:
   > - 3 counts from replicate 3, plus reaction, positive strand
   > - 1 count from replicate 4, minus reaction, negative strand
   > - 1 count from replicate 1, minus reaction, negative strand
   > - 1 count from replicate 3, plus reaction, negative strand
   - This result is for an experiment with mixed dinucleotides, for an experiment with a single dinucleotide, the average denominator is the total number of possible samples in the experiment that can give a signal (i.e. 4 for 2 replicates each of plus and minus reactions).
