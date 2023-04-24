# De novo assembly with SPAdes

# Map trimmed reads back to the assembled contigs using BWA-MEM

# This is a common practice in de novo assembly analysis to evaluate the assembly quality
# to perform variant calling (such as identifying SNPs). 
# Mapping the reads back to the contigs can provide useful information about the coverage, potential misassemblies, 
# and other metrics that can help assess the quality of the assembly.

echo "Mapping reads back to the contigs..."
BWA_EXEC="./bwa/bwa"
REF=assembly_output/contigs.fasta
$BWA_EXEC index $REF
$BWA_EXEC mem $REF trimmed_reads_1.fastq trimmed_reads_2.fastq > aligned_de_novo.sam

# Check the quality of the sam file 
samtools stats aligned_reads.sam > cross_assembly_stats.txt
echo "Assembly stats have been generated".

# The fact that a high proportion of reads (88,568 out of 97,602) are both properly paired and mapped is generally a good indication of the overall quality of the assembly.
# Properly paired reads are those for which the orientation and distance between the paired ends is consistent with the expected library preparation, which can help improve the accuracy and completeness of the assembly. 
# Mapped reads indicate that the sequencing reads were aligned to the assembly, which is a critical step in the assembly process.
# Additionally, the error rate is relatively low (0.27%) and the average quality score is high (37.5), both of which suggest good overall quality of the assembly.

# Now we can map it to a reference
# Have you copied the reference file to the current directory

######### A C T U A L    M A P P I N G   TO  R E F E R E N C E ##########

# intiating bwa from a directory
echo "Mapping reads back to the contigs..."
BWA_EXEC="./bwa/bwa"
REF=ref.fasta
$BWA_EXEC index $REF
$BWA_EXEC mem $REF ./assembly_output/contigs.fasta > aligned_to_ref.sam

# convert to binary
samtools view -Sb aligned_to_ref.sam | samtools sort -o aligned_to_ref.bam -

# index the sorted bam
samtools index aligned_to_ref.bam


# console_output 
# Mapping reads back to the contigs...
# [bwa_index] Pack FASTA... 0.00 sec
# [bwa_index] Construct BWT for the packed sequence...
# [bwa_index] 0.00 seconds elapse.
# [bwa_index] Update BWT... 0.00 sec
# [bwa_index] Pack forward-only FASTA... 0.00 sec
# [bwa_index] Construct SA from BWT and Occ... 0.01 sec
# [main] Version: 0.7.17-r1198-dirty
# [main] CMD: ./bwa/bwa index ref.fasta
# [main] Real time: 0.054 sec; CPU: 0.018 sec
# [M::bwa_idx_load_from_disk] read 0 ALT contigs
# [M::process] read 77 sequences (43211 bp)...
# [M::mem_process_seqs] Processed 77 reads in 0.035 CPU sec, 0.035 real sec
# [main] Version: 0.7.17-r1198-dirty
# [main] CMD: ./bwa/bwa mem ref.fasta ./assembly_output/contigs.fasta
# [main] Real time: 0.046 sec; CPU: 0.041 sec