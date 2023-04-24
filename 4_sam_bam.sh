#!/bin/bash

# Make sure that you have activate conda environment of samtools by using conda activate samtools_env

# This command converts a SAM file of aligned reads to a BAM file, which is a binary version of the same file format that is faster to process.
# After piping, sorts the aligned BAM file by the reference sequence position, which is required for downstream analysis such as variant calling and visualization.

samtools view -Sb aligned_to_ref.sam | samtools sort -o aligned_to_ref.bam -

# This command creates an index file for the sorted BAM file, which allows for faster random access to specific regions of the genome during analysis.
samtools index aligned_to_ref.bam

# Generate an mpileup file for varscan use
samtools mpileup -f ref.fasta aligned_to_ref.bam > sorted.mpileup

