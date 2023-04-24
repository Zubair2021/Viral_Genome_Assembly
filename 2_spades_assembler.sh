#!/bin/bash

# De novo assembly with SPAdes
echo "Assembling the genome..."

spades.py -t 4 -1 trimmed_reads_1.fastq -2 trimmed_reads_2.fastq -o assembly_output

echo "Check the directory assembly_output"


# This command runs the SPAdes genome assembler with specific input and output parameters.

# spades.py is the name of the SPAdes genome assembler executable file. By typing this command, you are instructing the command line to run this executable file.
# -t 4 specifies the number of threads that SPAdes should use for assembly. In this case, SPAdes will use four threads to run the assembly process.
# -1 trimmed_reads_1.fastq specifies the path to the first input file of paired-end reads. In this example, the file name is "trimmed_reads_1.fastq", and it contains the trimmed reads from the first end of a paired-end sequencing experiment. The -1 flag tells SPAdes that this is the first read file in a paired-end sequencing experiment.
# -2 trimmed_reads_2.fastq specifies the path to the second input file of paired-end reads. In this example, the file name is "trimmed_reads_2.fastq", and it contains the trimmed reads from the second end of a paired-end sequencing experiment. The -2 flag tells SPAdes that this is the second read file in a paired-end sequencing experiment.
# -o assembly_output specifies the output directory where the assembled genome will be saved. In this example, the output directory is named "assembly_output".


# Some of the options that could be used for spades.py
# --careful  # This option can be used to increase the accuracy of the assembly process by reducing the number of mismatches and indels in the resulting genome.

# --only-assembler  # This option can be used to skip the error correction and read filtering steps, which can speed up the assembly process for high-quality input data.

# --cov-cutoff auto  # This option can be used to specify a coverage cutoff for filtering low-coverage contigs from the final assembly.

# --meta  # This option can be used for metagenomic assembly, which requires additional options to handle multi-genome assembly.
