#!/bin/bash


# Varscan SNPs reported as vcf files
varscan mpileup2snp sorted.mpileup --min-reads2 10 --min-var-freq 0.01 --output-vcf 1 > variants.vcf


# mpileup2snp: This specifies that VarScan should identify SNPs from the mpileup output.
# sorted_al.bam: This specifies the name of the input BAM file containing the aligned sequencing reads.
# --min-reads2 10: This specifies that a minimum of 10 reads are required to call a SNP.
# --min-var-freq 0.01: This specifies that the minimum variant frequency (i.e., the proportion of reads supporting the variant) required to call a SNP is 0.01 (i.e., 1%).
# --output-vcf 1: This specifies that the output should be in the VCF format.
# > variants.vcf: This redirects the output to a file named "variants.vcf".