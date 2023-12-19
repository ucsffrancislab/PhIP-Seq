#!/bin/bash

#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH --mem=4G

module load gcc/6.2.0
module load bowtie/1.2.2
module load samtools/1.3.1

for fq in raw.data/*.fastq; do 

    bowtie -3 25 -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet path_to_vir3_reference_fasta_and_index_files/vir3 $fq | samtools view -u - | samtools sort -T ${fq%.fastq}.2.temp.bam -o ${fq%.fastq}.bam   
 
done

