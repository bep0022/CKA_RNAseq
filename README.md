# Cox Known Aging Anole RNAseq
## *This repository contains data and scripts for assessing RNAseq reads based on the combined effects of sex and age cohorts in Anolis sagrei.* 

## Data Accessibility 
### *"Primary paper associated with Data"*  
 **BioProject** ___ **DOI** [name for link](link)
### *rAnoSag1.mat Reference Genome and Annotation* 
  **Genome:** **BioProject** PRJNA1080665 **DOI** [name for link](link)
  
  **Annotation:** GCF_037176765.1-RS_2024_08 **BioProject** PRJNA1080665 **DOI** [name for link](link)

## Repository Outline and Summary
### Scripts
- 01_FastQC.sh - Use FastQC to evaluate the quality of the raw reads. Outputs "PreCleanQuality" folder/tarball of FastQC results.
  
- 02_Trim_FastQC.sh - Using .fastq sequencing reads downloaded, trim sequencing adapters and low quality regions of the reads with Trimmomatic; then use FastQC to evaluate the quality of the cleaned reads. Requires a file of adapter sequences used in sequencing. Outputs trimmed paired and unpaired reads and "PostCleanQuality" folder/tarball of FastQC results. 
  
- 03_mapping_HiSAT2.sh - Map the cleaned reads to the indexed reference genome, convert mapped reads .sam file to .bam using Samtools, use StringTie to create counts files of mapped reads, and convert counts files to counts matrices using python. Requires genome .fasta/.fna file and annotation .gff3/.gff file. 
   
- 04_.R - Use R script to create a .csv file of averaged normalized mapping counts across sexes and age cohorts in months to analyze target genes. 

### Output_Data
Directory containing output data from running the scripts above. Detailed descriptions are as follows.

- PreCleanQuality - Directory containing output from running script 1_Download_FastQC.sh. Contains zipped and .html links to FastQC quality assessment results for each SRR downloaded from NCBI. _1 refers to forward reads, _2 refers to reverse reads in the context of paired-end sequencing. 

- PostCleanQuality - Directory containing output from running script 2_Trim_FastQC.sh. Contains zipped and .html links to FastQC quality assessment results for each SRR downloaded from NCBI that were cleaned and trimmed. Like above, _1 refers to forward reads, _2 refers to reverse reads in the context of paired-end sequencing. paired refers to SRR's whose forward and reverse reads were kept after trimming. unpaired refers to those SRR's missing either the forward (_1) or reverse (_2).

- HiSAT2 - contains all input and output data that corresponds with the initial input of gene_count_matrix.csv created by 3a_mapping_HiSAT2.sh.
   
- STAR - contains all input and output data that corresponds with the gene_count_matrix.csv created by 3b_mapping_STAR.sh.

