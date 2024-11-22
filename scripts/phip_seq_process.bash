#!/usr/bin/env bash
#SBATCH --export=NONE   # required when using 'module'

hostname
echo "Slurm job id:${SLURM_JOBID}:"
date


set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail
if [ -n "$( declare -F module )" ] ; then
	echo "Loading required modules"
	module load CBI samtools r
fi
#set -x	#	print expanded command before executing it
set -v


#	initially based on /francislab/data1/working/20240925-Illumina-PhIP/20241009-PhIP/phip_seq_process.bash


function usage(){
	set +x
	echo
	echo "Usage:"
	echo
	echo $0 --manifest my_manifest.csv --output output/
	echo
	echo "Manifest CSV format:"
	echo "subject,sample,bam_file,input"
	echo
	exit
}



MANIFEST="/francislab/data1/raw/20240925-Illumina-PhIP/manifest.csv"
#	subject,sample,bam_file,input

OUTPUT="out"

while [ $# -gt 0 ] ; do
	case $1 in
		-m|--manifest)
			shift; MANIFEST=$1; shift;;
		-o|--output)
			shift; OUTPUT=$1; shift;;
		*)
			echo "Unknown param :${1}:"; usage ;;
	esac
done


mkdir -p ${OUTPUT}



echo "Counting"

mkdir -p ${OUTPUT}/counts

while read sample bam ; do
	echo $sample $bam
	if [ -f $bam ] ; then
		samtools view -F SECONDARY,SUPPLEMENTARY -q 40 ${bam} \
			| awk '( $4 < 10 ){print $3}' | sort -k1n | uniq -c \
			| awk '{print $2","$1}' | sed -e "1 i id,${sample}" \
			| gzip > ${OUTPUT}/counts/${sample}.q40.count.csv.gz
		chmod 440 ${OUTPUT}/counts/${sample}.q40.count.csv.gz
	fi
done < <( awk -F, '(NR>1){print $2"\t"$3}' ${MANIFEST} )




#	Note: To perform the Z-score analysis, count.combined files are merged into a table, and columns corresponding with no-serum controls are summed in a column called "input".

#	Separate the "input" samples

#	In this case, the dataset is 4 separate comparisono so do it separately

echo "Moving blanks into separate directory"

mkdir -p ${OUTPUT}/counts/input

while read sample ; do
	echo $sample
	if [ -f ${OUTPUT}/counts/${sample}.q40.count.csv.gz ] ; then
		mv ${OUTPUT}/counts/${sample}.q40.count.csv.gz ${OUTPUT}/counts/input
	fi
done < <( awk -F, '( $4 == "input" ){print $2}' ${MANIFEST} )







#	Then sum them with `sum_counts_files.py`


echo "Summing up blanks"

sum_counts_files.py -o ${OUTPUT}/input/All.count.csv ${OUTPUT}/input/*.q40.count.csv.gz
sed -i '1s/sum/input/' ${OUTPUT}/input/All.count.csv
gzip ${OUTPUT}/out/input/All.count.csv
chmod -w ${OUTPUT}/input/All.count.csv.gz






#	Create count matrices to feed to Z-score


echo "Merging in prep for Zscores"

merge_all_combined_counts_files.py --int -o ${OUTPUT}/All.count.csv ${OUTPUT}/*.q40.count.csv.gz ${OUTPUT}/input/All.count.csv.gz
chmod 400 ${OUTPUT}/All.count.csv






#	Full pipeline … Process as 3 replicates (a,b,c), 2 replicates ( (a,b),(a,c),(b,c))

#	Multiple runs with differing blank / input  and data for Zscore analysis - 4 different conditions


#	Really doesn't take that long so could just run locally.

#	Create Zscores

echo "Computing Zscores"

elledge_Zscore_analysis.R ${OUTPUT}/All.count.csv







#	Zscores with public epitopes





echo "Selecting public epitopes"


#	This depends on the number of samples.
#awk 'BEGIN{FS=OFS=","}{print $13,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${OUTPUT}/All.count.Zscores.csv > tmp
awk 'BEGIN{FS=OFS=","}{x=$NF-2;NF=NF-2;print x,$0}' ${OUTPUT}/All.count.Zscores.csv > tmp
#	L13,L14,L15,L19,L1,L20,L21,L2,L3,L7,L8,L9,id,group,input
#	id,L13,L14,L15,L19,L1,L20,L21,L2,L3,L7,L8,L9

head -1 tmp > ${OUTPUT}/All.count.Zscores.reordered.join_sorted.csv
tail -n +2 tmp | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.reordered.join_sorted.csv
join --header -t, /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.join_sorted.csv ${OUTPUT}/All.count.Zscores.reordered.join_sorted.csv > ${OUTPUT}/All.public_epitope_annotations.Zscores.csv 








#	I think that the output of the zscore command needs to have the columns reordered.

#	They are "all of the samples ...,id,group,input"

#	Not sure if "input" is needed anymore or what "group" is.

#	I don't use this.


#awk 'BEGIN{FS=OFS=","}(NR==1){print $5,$6,$1,$2,$3,$4,$7}(NR>1){printf "%d,%d,%.2g,%.2g,%.2g,%.2g,%d\n",$5,$6,$1,$2,$3,$4,$7}' \
#  Elledge/fastq_files/merged.combined.count.Zscores.csv > Elledge/fastq_files/merged.combined.count.Zscores.reordered.csv





#	Determine actual hits by zscore threshold in both replicates















echo "Booleanizing"

while read subject do
	echo $subject
#	if [ -f ${OUTPUT}/counts/${sample}.q40.count.csv.gz ] ; then
#		mv ${OUTPUT}/counts/${sample}.q40.count.csv.gz ${OUTPUT}/counts/input
#	fi

	all_samples=$( awk -F, -v subject=${subject} '( $1==subject ){print $2}' ${MANIFEST} | paste -sd' ' )
	echo $all_samples

	booleanize_Zscore_replicates.py --sample ${subject} --matrix ${OUTPUT}/All.count.Zscores.csv \
		--output ${OUTPUT}/${subject}.count.Zscores.hits.csv ${all_samples}

done < <( awk -F, '(NR>1){print $2}' ${MANIFEST} | sort | uniq )





#	**Account for public epitopes BEFORE virus score**


#	Create a list of viral species sorted by number of hits for all samples

echo "Creating a fixed species order"

tail -n +2 ${OUTPUT}/?.count.Zscores.hits.csv | sort -t, -k1,1 \
	| awk -F, '($2=="True")' > ${OUTPUT}/All.count.Zscores.merged_trues.csv

sed -i '1iid,all' ${OUTPUT}/All.count.Zscores.merged_trues.csv

join --header -t, ${OUTPUT}/All.count.Zscores.merged_trues.csv \
	/francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv > tmp

awk -F, '(NR>1){print $3}' tmp | sort | uniq -c | sort -k1nr,1 | sed 's/^ *//' \
	| cut -d' ' -f2- > ${OUTPUT}/species_order.txt
\rm tmp




#	Create virus scores in same order for all samples using the previous species order file.













echo "Calculate virus scores"

while read subject do
	echo $subject
	elledge_calc_scores_nofilter_forceorder.py --hits ${OUTPUT}/${subject}.count.Zscores.hits.csv \
		--oligo_metadata /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.csv \
		--species_order ${OUTPUT}/species_order.txt > tmp 
	head -1 tmp > ${OUTPUT}/${subject}.count.Zscores.hits.virus_scores.csv
	tail -n +2 tmp | sort -t, -k1,1 >> ${OUTPUT}/${subject}.count.Zscores.hits.virus_scores.csv
done < <( awk -F, '(NR>1){print $2}' ${MANIFEST} | sort | uniq )
merge_results.py --int -o ${OUTPUT}/merged.virus_scores.csv ${OUTPUT}/*.hits.virus_scores.csv








#	Threshold check


#	~/github/ucsffrancislab/PhIP-Seq/Elledge/VirScan_viral_thresholds.csv 


#	A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file “VirScan_viral_thresholds” contains the thresholds for each virus (Supplementary materials).


#	GREATER THAN THRESHOLD. NOT GREATER THAN OR EQUAL TO THE THRESHOLD.


echo "Thresholds"

for virus_scores in ${OUTPUT}/*.count.Zscores.hits.virus_scores.csv ; do
	echo ${virus_scores}
	join --header -t, ~/github/ucsffrancislab/PhIP-Seq/Elledge/VirScan_viral_thresholds.csv ${virus_scores} \
		| awk 'BEGIN{FS=OFS=","}(NR==1){print "Species",$3}(NR>1 && $3>$2){print $1,$3}' \
		> ${virus_scores%.csv}.threshold.csv
done 


















echo "Filter with join to public epitopes BEFORE"

for hits in ${OUTPUT}/*.count.Zscores.hits.csv ; do
	echo ${hits}
	join -t, <( tail -n +2 ${hits} | sort -t, -k1,1 ) \
		<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
		| awk -F, '($2=="True"){print $6}' | sort | uniq \
		> ${hits%.csv}.found_public_epitopes.BEFORE_scoring.txt
	sed -i '1iSpecies' ${hits%.csv}.found_public_epitopes.BEFORE_scoring.txt
done




echo "Filter with join to public epitopes AFTER"

#	*.7.peptides.txt files created during my version of the calc_virus_score scripts

for peptides in ${OUTPUT}/*.count.Zscores.hits.7.peptides.txt ; do
	echo ${peptides}
	join -t, <( sort -t, -k1,1 ${peptides} ) \
		<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
		| awk -F, '{print $2}' | sort | uniq \
		> ${peptides%.7.peptides.txt}.found_public_epitopes.AFTER_scoring.txt
	sed -i '1iSpecies' ${peptides%.7.peptides.txt}.found_public_epitopes.AFTER_scoring.txt
done








for scoring in ${OUTPUT}/*.count.Zscores.hits.found_public_epitopes.*_scoring.txt ; do
	echo ${scoring}
	join --header -t, ${scoring} ${scoring%.found_public_epitopes.*_scoring.txt}.virus_scores.threshold.csv \
		> ${scoring%.txt}.seropositive.csv
done 
merge_results.py --int -o ${OUTPUT}/merged.seropositive.csv ${OUTPUT}/*_scoring.seropositive.csv


