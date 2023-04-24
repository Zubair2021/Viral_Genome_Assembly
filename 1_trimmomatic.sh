#!/bin/bash

# Input files
READ1="S10_S10_L001_R1_001.fastq"
READ2="S10_S10_L001_R2_001.fastq"

# Trimming with Trimmomatic
echo "Trimming raw reads..."
TRIMMOMATIC_PATH="/home/zubair/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS_PATH="/home/zubair/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
java -jar $TRIMMOMATIC_PATH PE -phred33 $READ1 $READ2 trimmed_reads_1.fastq unpaired_1.fastq trimmed_reads_2.fastq unpaired_2.fastq ILLUMINACLIP:$ADAPTERS_PATH:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
