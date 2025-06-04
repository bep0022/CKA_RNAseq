# Cox Known Aging Anole RNAseq
## *This repository contains data and scripts for assessing RNAseq reads based on the combined effects of sex and age cohorts in Anolis sagrei.* 

## Data Accessibility 
### *"Towards complete and error-free genome assemblies of all vertebrate species"*  
 **BioProject** ___ **DOI** [ 10.1038/s41586-021-03451-0]( 10.1038/s41586-021-03451-0)
### *rAnoSag1.mat Reference Genome and Annotation: Completed by the Vertebrate Genome Project using samples obtained by Tonia S. Schwartz from Auburn University* 
  **Genome:** GCF_037176765.1-RS_2024_08
  **BioProject:** PRJNA1080665 **DOI:** [NCBI Link, for now](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_037176765.1/)
  
  **Annotation:** GCF_037176765.1-RS_2024_08 
  **BioProject:** PRJNA1080665 **DOI:** [NCBI Link, for now](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_037176765.1/)

## Repository Outline and Summary
### Scripts
- 01_FastQC.sh - Use FastQC to evaluate the quality of the raw reads. Outputs "PreCleanQuality" folder/tarball of FastQC results.
  
- 02_Trim_FastQC.sh - Using .fastq sequencing reads downloaded, trim sequencing adapters and low quality regions of the reads with Trimmomatic; then use FastQC to evaluate the quality of the cleaned reads. Requires a [file of adapter sequences](https://github.com/bep0022/CKA_RNAseq/blob/main/AdaptersToTrim_All.fa) used in sequencing. Outputs trimmed paired and unpaired reads and "PostCleanQuality" folder/tarball of FastQC results. 
  
- 03_Mapping_HiSAT2.sh - Map the cleaned reads to the indexed reference genome, convert mapped reads .sam file to .bam using Samtools, use StringTie to create counts files of mapped reads, and convert counts files to counts matrices using [python script](https://github.com/bep0022/CKA_RNAseq/blob/main/prepDE.py3). Requires genome .fasta/.fna file and annotation .gff3/.gff file. 
   
- 04_Normalize_Counts.R - Use R script to create a .csv file of averaged normalized mapping counts across sexes and age cohorts in months to analyze target genes. 

### Output_Data
Directory containing output data from running the scripts above. Detailed descriptions are as follows.

- PreCleanQuality - Directory containing output from running script 1_FastQC.sh. Contains zipped and .html links to FastQC quality assessment results for each sample. _1 refers to forward reads, _2 refers to reverse reads in the context of paired-end sequencing. 

- PostCleanQuality - Directory containing output from running script 2_Trim_FastQC.sh. Contains zipped and .html links to FastQC quality assessment results for each sample that were cleaned and trimmed. Like above, _1 refers to forward reads, _2 refers to reverse reads in the context of paired-end sequencing. paired refers to samples whose forward and reverse reads were kept after trimming. unpaired refers to those samples missing either the forward (_1) or reverse (_2).

- HiSAT2 - contains all output data that corresponds with the initial input of [gene_count_matrix.csv](https://github.com/bep0022/CKA_RNAseq/blob/main/Output_Data/HiSAT2/gene_count_matrix.csv) created by 3_Mapping_HiSAT2.sh, including the [avg_normalized_gene_count_matrix.csv](https://github.com/bep0022/CKA_RNAseq/blob/main/Output_Data/HiSAT2/avg_normalized_gene_count_matrix.csv) created by 04_Normalize_Counts.R.
   
