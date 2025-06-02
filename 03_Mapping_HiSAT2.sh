#!/bin/sh

#SBATCH --job-name=Mapping3_HiSat2                   # job name
#SBATCH --nodes=2                                    # node(s) required for job
#SBATCH --ntasks=10                                  # number of tasks across all nodes
#SBATCH --mem=120G                                   # amount of memory
#SBATCH --partition=general                          # name of partition
#SBATCH --time=48:00:00                              # Run time (D-HH:MM:SS)
#SBATCH --output=Maptest3-%j.out                     # Output file. %j is replaced with job ID
#SBATCH --error=Maptest3_error-%j.err                # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu               # CHANGE THIS EMAIL ADDRESS WITH YOUR EMAIL ADDRESS

##This script employs Easley through Auburn University; to run this script, use "sbatch [script]"
##see https://hpc.auburn.edu/hpc/docs/hpcdocs/build/html/easley/easley.html for more information

######## Script 3: Map cleaned .fastq reads to a reference genome using HISAT2 mapper ############
## Purpose: The purpose of this script is to map the cleaned RNASeq data (.fastq) to the reference genome.

## The reference genome for our study is rAnoSag1.mat (NCBI BioProject Accession PRJNA1080665) curated by The Vertebrate Genome Project,
## with annotation also created by The Vertebrate Genome Project (NCBI RefSeqGCF_037176765.1).
## The animals sequenced for the reference came from Tonia S. Schwartz's colony at Auburn University
## Genome Citation:

## Use HISAT2 to index the reference genome and map cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft format
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              Use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HISAT2 Indexing
##  Input: Reference genome file (.fna), and annotation file (.gff)
##  Output: Indexed genome
## HISAT2 Mapping
##  Input: Cleaned read files, paired (.fq.gz) created during step 2 (trimmomatics and FASTQC script); Indexed genome
##  Output: Alignment .sam files
## Samtools .sam to .bam conversion
##  Input: Alignment files .sam
##  Output: Sorted alignment .bam files
## Stringtie Counting reads
##  Input: sorted alignment .bam files
##  Output: Directories of counts files for Ballgown (R program for DGE)
## prepDE.py Python script to create counts matrices from the Stringtie output.
##  Input: Directories from Stringtie
##  Output:  .csv files of counts matrix
###############################################

module load hisat2/2.2.1
module load stringtie/2.1.6
module load python/2.7.1
module load gcc/9.3.0
module load samtools/1.19
module load bcftools/1.17
module load gffreader/12.7

#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories

WD=/scratch/<path>/<to>/<directory/RNAseq
CLEAND=/scratch/bep0022/01.RawData/CleanData                                                  ##This is where the cleaned, paired files are located
REFD=/home/bep0022/ReferenceGenomes/rAnoSag1.mat/ncbi_dataset/data/GCF_037176765.1            ##This is where the indexed reference genome is located
MAPD=/scratch/bep0022/01.RawData/Map_HiSat2                                                   ##This is where the mapped sequence files will be located (.sam and .bam)
COUNTSD=/scratch/bep0022/01.RawData/Counts_HiSat2_StringTie                                   ##This is where the read counts files will be located alongside other stringtie outputs

RESULTSD=/home/bep0022/CKA_RNAseq/Counts_HiSat2_Stringtie_Results                             ##This is where the results you want to be brought back to your computer are

REF=GCF_037176765.1_rAnoSag1.mat_genomic                                                      ## This is the shortcut name for the genome reference

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD
cp /home/bep0022/ReferenceGenomes/rAnoSag1.mat/ncbi_dataset/data/GCF_037176765.1/$REF.fna .
cp /home/bep0022/ReferenceGenomes/rAnoSag1.mat/ncbi_dataset/data/GCF_037176765.1/genomic.gff .

###  Identify exons and splice sites
gffread genomic.gff -T -o ${REF}.gtf               ## gffread -T -o converts the annotation file from .gff to .gft format for HiSat2 to use.
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build its own index.
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Asag_mat_index
#hisat2-build -p 6  ${REF}.fna Asag_mat_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd $CLEAND  #### This is where our cleaned paired reads are located.

## Create list of fastq files to map.
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fq.gz" |cut -d "_" -f 1,2| sort | uniq > list

## Move to the directory for mapping
cd $MAPD

## move the list of unique IDs from the original files to map
mv $CLEAND/list .

while read i;
do
  ## HiSat2 is the mapping program
  ## -p indicates number of processors
  ## --dta reports alignments for StringTie
  ## --rf is the read orientation
  ## -x basename of the index for the reference genome
  ## -1 Comma-separated list of files containing mate 1s
  ## -2 Comma-separated list of files containing mate 2s
  ## -S File to write SAM alignments to

hisat2 -p 6 --dta --phred33	  \
    -x "$REFD"/Asag_mat_index       \
    -1 "$CLEAND"/"$i"_1_paired.fq.gz  -2 "$CLEAND"/"$i"_2_paired.fq.gz      \
    -S "$i".sam

    ### view convert the SAM file into a BAM file
    ### -bS BAM is the binary format corresponding to the SAM text format.
    ### sort convert the BAM file to a sorted BAM file.
    ### flagstat counts the number of alignments for each FLAG type
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam

  samtools view -@ 6 -bS "$i".sam > "$i".bam

  samtools sort -@ 6 "$i".bam -o "$i"_sorted.bam


 ### Index the BAM and get mapping statistics
  samtools flagstat "$i"_sorted.bam > "$i"_Stats.txt



  ## Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model.
  ## -p Specify the number of processing threads (CPUs) to use for transcript assembly
  ## -e limits the processing of read alignments to estimating the coverage of the transcripts given with the -G option
  ## -B enables the output of Ballgown input table files
  ## -G Use a reference annotation file (in GTF or GFF3 format) to guide the assembly process
  ## -o Sets the name of the output GTF file where StringTie will write the assembled transcripts
  ## -l Sets "$i" as the prefix for the name of the output transcripts


mkdir "${COUNTSD}"/"$i"
stringtie -p 6 -e -B -G "$REFD"/"$REF".gtf -o "$COUNTSD"/"$i"/"$i".gtf -l "$i"   "$MAPD"/"$i"_sorted.bam
done<list

#####################  Copy Results from your mapping directory to HOME Directory. These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt $RESULTSD

### The prepDE.py3 is a python script that converts the files in your ballgown folder to a count matrix
cd $COUNTSD
python /home/bep0022/prepDE.py3 -i $COUNTSD

### copy the final results files (the count matricies that are .csv to your home directory)
cp *.csv $RESULTSD

#Output gene_count_matrix.csv from this script can be found in Output_Data/HiSat2
