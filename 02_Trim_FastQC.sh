#! /bin/bash

#SBATCH --job-name=Test3_Trim_FastQC                 # job name
#SBATCH --nodes=2                                    # node(s) required for job
#SBATCH --ntasks=10                                   # number of tasks across all nodes
#SBATCH --partition=general                          # name of partition
#SBATCH --time=24:00:00                              # Run time (D-HH:MM:SS)
#SBATCH --output=test3trim-%j.out                        # Output file. %j is replaced with job ID
#SBATCH --error=test3trim_error-%j.err                   # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu               # CHANGE THIS EMAIL ADDRESS WITH YOUR EMAIL ADDRESS


##This script employs Easley through Auburn University; to run this script, use "sbatch [script]"
##see https://hpc.auburn.edu/hpc/docs/hpcdocs/build/html/easley/easley.html for more information

######## Script 2: Trim SRA sequences as needed and perform FastQC ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##	 and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## FASTQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
##              Input Data: trimmed R1 and R2 .fastq reads
##              Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics. See Output_Data/PreCleanQuality for example.
##                      The last line of this script will make a tarball of the output directory to bring back to your computer
###############################################

# Modules
## to see current modules available, cd /tools
module load trimmomatic/0.39
module load fastqc/0.11.6

##########  Define variables and make directories
## Replace variable pathways with your specific information


# Variables: working directory(WD), raw data directory (DD), cleaned data directory (CD), post-cleaned status directory (CS), name of file containing the adpaters (adapters).
WD=/scratch/bep0022/01.RawData                                  ## Contains Raw and Cleaned (where we perform trimmomatics) Data directories
DD=/scratch/bep0022/01.RawData/All_Raw                          ## Contains raw fastq reads downloaded from NCBI
CD=/scratch/bep0022/01.RawData/CleanData                        ## Contains trimmed and cleaned fastq reads from running Trimmomatic
CS=PostCleanQuality
adapters=AdaptersToTrim_All.fa          ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. Can be found in main GitHub repository.
                                        ## You will likely need to edit this for other projects based on how your libraries
                                        ## were made to search for the correct adapters for your project.

## make the directories to hold the Cleaned Data (CD) files, and the directory to hold the results for assessing quality of the cleaned data (CS); note that your WD and DD
## should have already been made while running script 1_Download_FastQC.sh
mkdir $CD
mkdir $WD/$CS

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd $DD          #contains raw .fastq files from NCBI

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list,
        ## grep means grab all the file names that end in ".fastq",
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fq.gz" |cut -d "_" -f 1,2 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters).
        ## CHECK: You will need to edit this path for the file that is in the directory of interest containing said adapters file.
        ## Download Adapters file from GitHub using the command:
        ## scp [pwd]/AdaptersToTrim_All.fa [your ID]@easley.auburn.edu:[pwd to your desired easley directory]

cp /home/bep0022/AdaptersToTrim_All.fa .

### Run a while loop to process through the names in the list and trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format.
        ### Make sure to change pathway to access trimmomatic if not running on auburn easley.

        java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 6 -phred33 \
        "$i"_1.fq.gz "$i"_2.fq.gz  \
        $CD/"$i"_1_paired.fq.gz $CD/"$i"_1_unpaired.fq.gz  $CD/"$i"_2_paired.fq.gz $CD/"$i"_2_unpaired.fq.gz \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
                ## requiredQuality: specifies the average quality required.

        ############## FASTQC to assess quality of the Cleaned sequence data
        ## FastQC: run on each of the data files that have 'All' to check the quality of the data
        ## The output from this analysis is a folder of results and a zipped file of results

fastqc $CD/"$i"_1_paired.fq.gz --outdir=$WD/$CS
fastqc $CD/"$i"_2_paired.fq.gz --outdir=$WD/$CS

done<list                       # This is the end of the loop

## move to the directory with the cleaned data
cd $WD/$CS

####### Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate
tar cvzf $CS.tar.gz $WD/$CS/*

# All output files from this script can be found in Output_Data/PostCleanQuality
