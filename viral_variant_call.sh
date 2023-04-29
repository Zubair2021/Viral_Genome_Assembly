#!/bin/bash

##############################################################################################
#                   V A R I A N T    C A L L I N G    P I P E L I N E 
# --------------------------------------------------------------------------------------------
# Script for variant calling using BWA-MEM and SAMtools on the Alabama Supercomputer (ASC).     #
# Author: Zubair Khalid (zubair.khalid@auburn.edu)
# Dated:  04/29/22
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                   #
# Input:                                                                                     #
#     1. Directory containing the raw sequencing files                                        #
#     2. Path to the reference genome file (in fasta format)                                  #
#                                                                                            #
# Output:                                                                                    #
#     1. final_bam_files - Directory containing the final sorted BAM files                   #
#     2. contigs_files - Directory containing the assembled contigs                          #
#     3. assembly_stats_files - Directory containing assembly statistics files                #
#     4. sam_files - Directory containing the SAM files                                      #
#     5. mpileup_files - Directory containing the mpileup files for Varscan analysis          #
#     6. bcfs - Directory containing the BCF files for each sample                            #
#                                                                                            #
# Usage:                                                                                     #
#     bash variant_calling.sh input_dir ref.fasta                                             #
#                                                                                            #
# Note:                                                                                      #
#     - Ensure that all required modules are loaded (trimmomatic, spades, bwa, samtools)     #
#     - Ensure that ref.fasta is present in the current directory                             #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################################################

# While running the script, specify the directory containing the raw fastq seqs as command argument 1
echo "Please enter the name/path of the raw fastq directory:"
read -e INPUT_DIR

# Tell user if the file is missing
if [ ! -d $INPUT_DIR ]; then
    echo "Error: Please enter the corrent name/path of the directory containing fastq files"
    exit 1
fi

# Check if the reference.fasta has been added to the directory
# Prompt the user to ensure ref.fasta is in the current directory
echo "Please enter the name of fasta file:"
read REF

# Tell user if the file is missing
if [ ! -f $REF ]; then
    echo "Error: reference fasta file is missing from the current directory."
    echo "Please copy reference fasta to the current directory and run the script again."
    exit 1
fi

# Load required modules (these are installed on ASC by default)
# One can check the availability & path of a module by using module avail <toolname>
module load trimmomatic
module load spades
module load bwa
module load samtools
module load bcftools

# Set the output directory for the final sorted BAM files
FINAL_BAM_DIR="final_bam_files"
CONTIGS_DIR="contigs_files"
ASSEMBLY_STATS_DIR="assembly_stats_files"
SAM_DIR="sam_files"
VCF_DIR="vcf_files"

# Create the final BAM directory if it doesn't exist
if [ ! -d $FINAL_BAM_DIR ]; then
    mkdir $FINAL_BAM_DIR
fi

# Create the contigs directory if it doesn't exist
if [ ! -d $CONTIGS_DIR ]; then
    mkdir $CONTIGS_DIR
fi

# Create the assembly stats directory if it doesn't exist
if [ ! -d $ASSEMBLY_STATS_DIR ]; then
    mkdir $ASSEMBLY_STATS_DIR
fi

# Create the SAM directory if it doesn't exist
if [ ! -d $SAM_DIR ]; then
    mkdir $SAM_DIR
fi

# Create the VCF directory if it doesn't exist
if [ ! -d $VCF_DIR ]; then
    mkdir $VCF_DIR
fi

# Loop over all read pairs in the input directory and run Trimmomatic, SPAdes, and BWA-MEM
for READ1 in $INPUT_DIR/*_R1_*.fastq
do
    # Check if the file exists and skip if not
    if [ ! -e $READ1 ]; then
        continue
    fi

    # Extract the sample name from the READ1 filename
    SAMPLE=$(basename $READ1 | cut -d "_" -f 1)

    # Define input and output filenames
    READ2=${READ1/_R1_/_R2_}
    TRIMMED_READ1=${SAMPLE}_trimmed_reads_1.fastq
    TRIMMED_READ2=${SAMPLE}_trimmed_reads_2.fastq
    OUTPUT_DIR=${SAMPLE}_assembly_output
    CONTIGS=${SAMPLE}_contigs.fasta
    SCAFFOLDS=${SAMPLE}_scaffolds.fasta
    ALIGNED_DE_NOVO=${SAMPLE}_aligned_de_novo.sam
    CROSSLINK_STATS=${SAMPLE}_cross_assembly_stats.txt
    ALIGNED_TO_REF=${SAMPLE}_aligned_to_ref.sam
    SORTED_BAM=${SAMPLE}_aligned_to_ref.sorted.bam
    MPILEUP=${SAMPLE}_sorted.mpileup
    VCF=${SAMPLE}.vcf

    # Create the output directory if it doesn't exist
    if [ ! -d $OUTPUT_DIR ]; then
         mkdir $OUTPUT_DIR
    fi

    # Run Trimmomatic with the given parameters and output to the output directory
    echo "Trimming $READ1 and $READ2..."
    trimmomatic PE -phred33 $READ1 $READ2 $OUTPUT_DIR/$TRIMMED_READ1 $OUTPUT_DIR/unpaired_1.fastq $OUTPUT_DIR/$TRIMMED_READ2 $OUTPUT_DIR/unpaired_2.fastq TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run SPAdes with the given parameters and output to the output directory
    echo "Assembling $TRIMMED_READ1 and $TRIMMED_READ2..."
    spades.py -t 4 -1 $OUTPUT_DIR/$TRIMMED_READ1 -2 $OUTPUT_DIR/$TRIMMED_READ2 -o $OUTPUT_DIR

    mv $OUTPUT_DIR/contigs.fasta $OUTPUT_DIR/$CONTIGS
    mv $OUTPUT_DIR/scaffolds.fasta $OUTPUT_DIR/$SCAFFOLDS

    # Check the output directory
    echo "Check the directory $OUTPUT_DIR"

    # Map trimmed reads back to the assembled contigs using BWA-MEM
    echo "Mapping reads back to the contigs..."
    BWA_EXEC="bwa"
    
    REF_CONTIGS=$OUTPUT_DIR/$CONTIGS
    $BWA_EXEC index $REF_CONTIGS
    $BWA_EXEC mem $REF_CONTIGS $OUTPUT_DIR/$TRIMMED_READ1 $OUTPUT_DIR/$TRIMMED_READ2 > $OUTPUT_DIR/$ALIGNED_DE_NOVO

    # Check the quality of the sam file 
    samtools stats $OUTPUT_DIR/$ALIGNED_DE_NOVO > $OUTPUT_DIR/$CROSSLINK_STATS
    echo "Assembly stats have been generated."

    # Map reads back to a reference genome
    echo "Mapping reads back to the reference genome..."
    $BWA_EXEC index $REF
    $BWA_EXEC mem $REF $OUTPUT_DIR/$CONTIGS > $OUTPUT_DIR/$ALIGNED_TO_REF
    

    # Convert the sam file to bam format
    echo "Converting to BAM format..."
    samtools view -Sb $OUTPUT_DIR/$ALIGNED_TO_REF > $OUTPUT_DIR/$SORTED_BAM

    # Sort the BAM file and create an index
    echo "Sorting the BAM file and creating an index..."
    samtools sort -o $FINAL_BAM_DIR/$SORTED_BAM $OUTPUT_DIR/$SORTED_BAM
    
    # To check the order of the genes/chr in the sorted bam
    # samtools idxstats S11_aligned_to_ref.sorted.bam | cut -f 1 | sort | uniq
    # In the console, each row represents one uniq gene, ideally the reference should be ordered alphabetically as well

    samtools index $FINAL_BAM_DIR/$SORTED_BAM

    # Generate an mpileup file for Varscan
    echo "Generating mpileup file for Varscan..."
    samtools mpileup -f $REF $FINAL_BAM_DIR/$SORTED_BAM > $FINAL_BAM_DIR/$MPILEUP
    
    # Call variants with bcftools
    echo "Calling variants with bcftools..."
 
    bcftools mpileup -Ou -f $REF $FINAL_BAM_DIR/$SORTED_BAM | bcftools call -mv -Ov -o $VCF_DIR/$VCF

    # Move the generated files to their respective directories
    mv $OUTPUT_DIR/$CONTIGS $CONTIGS_DIR/$CONTIGS
    mv $OUTPUT_DIR/$CROSSLINK_STATS $ASSEMBLY_STATS_DIR/$CROSSLINK_STATS
    mv $OUTPUT_DIR/$ALIGNED_DE_NOVO $SAM_DIR/$ALIGNED_DE_NOVO
    mv $OUTPUT_DIR/$ALIGNED_TO_REF $SAM_DIR/$ALIGNED_TO_REF

done

ls -l vcf_files