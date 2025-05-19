#! /bin/bash
#SBATCH --job-name=Test5_CKA_FastQC            # job name
#SBATCH --nodes=1                                    # node(s) required for job
#SBATCH --ntasks=1                                   # number of tasks across all nodes
#SBATCH --partition=general                          # name of partition
#SBATCH --time=24:00:00                              # Run time (D-HH:MM:SS)
#SBATCH --output=test5-%j.out                        # Output file. %j is replaced with job ID
#SBATCH --error=test5_error-%j.err                   # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu               # CHANGE THIS EMAIL ADDRESS WITH YOUR EMAIL ADDRESS

##This script employs Easley through Auburn University; to run this script, use "sbatch [script]"
##see https://hpc.auburn.edu/hpc/docs/hpcdocs/build/html/easley/easley.html for more information

######### Script 1: Perform Preliminary FastQC ############
## Purpose:
##	1. use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

## FASTQC Pipeline:
##  Input: Downloaded SRA files .fastq
##  Output: Folder for each input .fastq SRA file; The last line of this script will make a tarball of the output directory (PreCleanQuality) to bring back to your computer

###############################################


########## Load Modules
module load fastqc/0.11.6

##########  Define variables and make directories
## Replace variable pathways with your specific information


## Make variables that represent your Raw data directory (DD) and your working directory (WD) in scratch, and the precleaned status directory (CS) to be tar-balled at the end.
DD=/scratch/bep0022/01.RawData/All_Raw                 #this directory is where your fastq files downloaded from NCBI will be dumped using fastq-dump
CS=/scratch/bep0022/01.RawData/PreCleanQuality         #this directory lies within your working directory and will contain the FastQC output files

#mkdir -p $DD
##  make the directories in SCRATCH for holding the raw data
##  -p tells it to make any upper level directories that are not there.
mkdir $CS

cd $DD

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results

fastqc *.fq.gz --outdir=$WD/$CS

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate (your CS directory)
cd $WD/$CS
tar -cvzf $CS.tar.gz $WD/$CS/*

#####Bringing tarball to home directory
#using the command scp [easleyID]@easley.auburn.edu:[pwd to tarball]/PreCleanQuality.tar.gz [directory to move the tarball to]

#All output data for this project can be found in Output_Data/PreCleanQuality
