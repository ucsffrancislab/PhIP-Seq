#!/usr/bin/env bash

#	From section "Data analysis" A of "VirScan: High-throughput Profiling of Antiviral Antibody Epitopes"

set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail
if [ -n "$( declare -F module )" ] ; then
	echo "Loading required modules"
	#module load CBI samtools
	module load bowtie	#/1.2.2
	module load samtools	#/1.3.1
fi
set -x


SAMPLE_ID=SAMPLE_ID
INDEX=~/github/ucsffrancislab/PhIP-Seq/Elledge/vir3

while [ $# -gt 1 ] ; do
	case $1 in
		-x|-i|--index)
			shift; INDEX=$1; shift;;
		-s|--sample_id)
			shift; SAMPLE_ID=$1; shift;;
		*)
			echo "Unknown param :${1}:"; exit 1 ;;
	esac
done

fq=$1





# aligning with bowtie like this effective requires a near exact, end-to-end forward match.
#	bowtie2 could do this as well. not sure if it would be better or faster. just different.






#	Not sure why the `samtools view -u` is needed.
#	  -u, --uncompressed         Uncompressed BAM output (and default to --bam)

bowtie -3 25 -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet \
	-x ${INDEX} \
	$fq \
	| samtools view -u - \
	| samtools sort -T ${fq%.fastq.gz}.2.temp.bam -o ${fq%.fastq.gz}.bam


#	  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
#	  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
#	  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
#	  -l/--seedlen <int> seed length for -n (default: 28)
#	  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
#	  --nofw/--norc      do not align to forward/reverse-complement reference strand
#	  -y/--tryhard       try hard to find valid alignments, at the expense of speed
#	  --best             hits guaranteed best stratum; ties broken by quality
#	  --quiet            print nothing but the alignments
#	  -S/--sam           write hits in SAM format
#	  -p/--threads <int> number of alignment threads to launch (default: 1)


#3. Check the alignment report file that ends in ".out"
#Note: Typically, >85% of the reads align to the reference file.


bam=${fq%.fastq.gz}.bam
samtools index $bam


#5. Count indexes with the following commands. The output is a file that ends in ".count.csv"



#	perhaps add SAMPLE_ID as an incoming parameter



#samtools idxstats $bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > ${bam%.bam}.count.csv
samtools idxstats $bam | cut -f 1,3 | sed -e '/^\*\t/d' -e "1 i id\t${SAMPLE_ID}" | tr "\\t" "," > ${bam%.bam}.count.csv


csv=${bam%.bam}.count.csv
gzip $csv



#	Why? Doesn't seem to be used.
#mkdir log_directory



#	8. If the same sample is run on two or more lanes of a flow cell and separate files are provided for each flow cell, combine the counts files from the different lanes using the following commands. These commands require the python script "combine_two_lanes.py" to be copied to the folder where you are running the commands (Supplementary materials).
#	
#	Note: In the code below, the samples were run on four lanes of an Illumina Nextseq 500 flow cell. The suffix of each count file is "L001_R1_001.count.csv.gz" if the count file was from the first lane of the flow cell, "L002_R1_001.count.csv.gz" if the count was from the second lane of the flow cell, etc.
#	
#	```
#	module load gcc/6.2.0
#	
#	module load python
#	
#	for i in raw.data/*L001_R1_001.count.csv.gz; do python combine_two_lanes.py $i
#	${i%1_R1_001.count.csv.gz}2_R1_001.count.csv.gz 
#	${i%1_R1_001.count.csv.gz}1_2_R1_001.count.csv; done
#	
#	for i in raw.data/*L003_R1_001.count.csv.gz; do python combine_two_lanes.py $i
#	${i%3_R1_001.count.csv.gz}4_R1_001.count.csv.gz 
#	${i%3_R1_001.count.csv.gz}3_4_R1_001.count.csv; done
#	
#	for i in raw.data/*L001_2_R1_001.count.csv; do python combine_two_lanes.py $i
#	${i%1_2_R1_001.count.csv}3_4_R1_001.count.csv 
#	${i%1_2_R1_001.count.csv}1_2_3_4_R1_001.count.combined.csv; done
#	```






#	9. Gzip the count.combined files with the following command.
#	
#	```
#	for i in raw.data/*1_2_3_4_R1.count.combined.csv; do gzip $i; done
#	```





