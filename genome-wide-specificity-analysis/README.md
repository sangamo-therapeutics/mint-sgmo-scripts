# README for processing all genome-wide specificity analysis code

The code here is required to generate all genome-wide specificity data in the paper (all data associated with Supplementary Tables 6, 8-13, 23-29, and 30-36).

There are 4 iterations of this code to account for experimental differences. All account for files with the following properties:

> - Each biological sample is analyzed using separate "plus" and "minus" primers in the assay that correspond to either direction of amplification with respect to the donor or on-target genomic integration site (depending on the assay). In these cases, "plus" refers to upstream of the dinucleotide and "minus" downstream for either the donor, for the genome-wide integration site nomination assay, or on-target genomic integration site for the genome-wide translocation nomination assay.
> - If using replicates, each file has "rep#" in it's name to correspond to which one.
> - Files required for the analysis are provided in separate R1 and R2 format in the SRA and are required to be in unpaired format to run the analysis scripts.
> - Make sure all files associated with a set of experimental conditions (same reagents and dosing) are all named the same except for these variable names (plus/minus and rep#) and do not have any other instance of these strings in the filenames. For example, here are all necessary files associated with Supplementary Table 10:
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_plus_rep1_R1.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_plus_rep1_R2.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_plus_rep2_R1.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_plus_rep2_R2.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_minus_rep1_R1.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_minus_rep1_R2.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_minus_rep2_R1.fastq.gz
>   - SI_TABLE_10_MACO1L_MACO1R_bxb1_minus_rep2_R2.fastq.gz

`pipeline_mixednt_5_reps.sh` and `pipeline_mixednt_1_rep.sh` both are used for datasets that were generated using mixed dinucleotide donors.
- `pipeline_mixednt_5_reps.sh` is used on data associated with Supplementary Table 6 and corresponds to processing and summarizing candidate integration sites across 5 replicates
- `pipeline_mixednt_1_reps.sh` is used on data associated with Supplementary Tables 8-9 and corresponds to processing and summarizing candidate integration sites with only 1 replicate.

`pipeline_nonpalindromic_nonmixed_2_reps.sh` is used datasets that used a single dinucleotide donor and 2 replicates.

The docker container set-up should be performed once and will allow the execution of all specificity scripts.

- Docker instructions:

