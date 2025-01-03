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
	echo "subject,sample,bam_file,type"
	echo
	exit
}



#MANIFEST="/francislab/data1/raw/20240925-Illumina-PhIP/manifest.csv"
#	subject,sample,bam_file,input

OUTPUT="out"
Q=40
SPECIES_ORDER=""
STOP_AFTER_ZSCORE=false
THRESHOLDS=3.5

while [ $# -gt 0 ] ; do
	case $1 in
		-m|--manifest)
			shift; MANIFEST=$1; shift;;
		-o|--output)
			shift; OUTPUT=$1; shift;;
		-s|--species|--species_order)
			shift; SPECIES_ORDER=$1; shift;;
		-t|--threshold*)
			shift; THRESHOLDS=$( echo $1 | tr , " " ); shift;;
		--stop_after_zscore)
			STOP_AFTER_ZSCORE=true; shift;;
		-q)
			shift; Q=$1; shift;;
		*)
			echo "Unknown param :${1}:"; usage ;;
	esac
done


mkdir -p ${OUTPUT}



echo "Counting"

mkdir -p ${OUTPUT}/counts/input

while read sample bam type ; do
	echo ":${sample}:${bam}:${type}:"

	if [ "${type}" == "input" ] ; then
	#	echo "true"
		f=${OUTPUT}/counts/input/${sample}.q${Q}.count.csv.gz
	else
	#	echo "false"
		f=${OUTPUT}/counts/${sample}.q${Q}.count.csv.gz
	fi

	#echo "after"

	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		if [ -f $bam ] ; then
			samtools view -F UNMAP,SECONDARY,SUPPLEMENTARY -q ${Q} ${bam} \
				| awk '( $4 < 10 ){print $3}' | sort -k1n | uniq -c \
				| awk '{print $2","$1}' | sed -e "1 i id,${sample}" \
				| gzip > ${f}
			chmod 440 ${f}
		fi
	fi
done < <( awk -F, '(NR>1){print $2"\t"$3"\t"$4}' ${MANIFEST} )




#	Note: To perform the Z-score analysis, count.combined files are merged into a table, and columns corresponding with no-serum controls are summed in a column called "input".

#	Separate the "input" samples

#	In this case, the dataset is 4 separate comparisono so do it separately
#
#	echo "Moving blanks into separate directory"
#
#	mkdir -p ${OUTPUT}/counts/input
#
#	while read sample ; do
#		echo $sample
#		f=${OUTPUT}/counts/input/${sample}.q40.count.csv.gz
#		if [ -f ${f} ] && [ ! -w ${f} ] ; then
#			echo "Write-protected ${f} exists. Skipping."
#		else
#	#		if [ -f ${OUTPUT}/counts/${sample}.q40.count.csv.gz ] ; then
#			mv ${OUTPUT}/counts/${sample}.q40.count.csv.gz ${OUTPUT}/counts/input
#			chmod -w ${f}
#	#		fi
#		fi
#	done < <( awk -F, '( $4 == "input" ){print $2}' ${MANIFEST} )







#	Then sum them with `sum_counts_files.py`


echo "Summing up blanks"

f=${OUTPUT}/counts/input/All.count.csv.gz
if [ -f ${f} ] && [ ! -w ${f} ] ; then
	echo "Write-protected ${f} exists. Skipping."
else
	sum_counts_files.py --int -o ${OUTPUT}/counts/input/All.count.csv ${OUTPUT}/counts/input/*.q${Q}.count.csv.gz
	sed -i '1s/sum/input/' ${OUTPUT}/counts/input/All.count.csv
	gzip ${OUTPUT}/counts/input/All.count.csv
	chmod -w ${f}
fi






#	Create count matrices to feed to Z-score


echo "Merging in prep for Zscores"

f=${OUTPUT}/All.count.csv
if [ -f ${f} ] && [ ! -w ${f} ] ; then
	echo "Write-protected ${f} exists. Skipping."
else
	merge_all_combined_counts_files.py --de_nan --int -o ${f} \
		${OUTPUT}/counts/*.q${Q}.count.csv.gz ${OUTPUT}/counts/input/All.count.csv.gz
	chmod -w ${f}
fi






#	Full pipeline … Process as 3 replicates (a,b,c), 2 replicates ( (a,b),(a,c),(b,c))

#	Multiple runs with differing blank / input  and data for Zscore analysis - 4 different conditions


#	Really doesn't take that long so could just run locally.

#	Create Zscores

echo "Computing Zscores"

f=${OUTPUT}/All.count.Zscores.csv
if [ -f ${f} ] && [ ! -w ${f} ] ; then
	echo "Write-protected ${f} exists. Skipping."
else

	#	elledge_Zscore_analysis.R ${OUTPUT}/All.count.csv
	#	mv ${OUTPUT}/All.count.Zscores.csv ${OUTPUT}/tmp1.csv
	#	awk 'BEGIN{FS=OFS=","}{x=$(NF-2);NF=NF-3;print x,$0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv


	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp1.csv #--count_threshold 200
	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv
	mv ${OUTPUT}/tmp1.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.csv

	head -1 ${OUTPUT}/tmp2.csv > ${OUTPUT}/All.count.Zscores.csv
	tail -n +2 ${OUTPUT}/tmp2.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.csv

	chmod -w ${f}
	#\rm ${OUTPUT}/tmp1.csv ${OUTPUT}/tmp2.csv
	\rm ${OUTPUT}/tmp2.csv

fi






#	Zscores with public epitopes

echo "Selecting public epitopes"
f=${OUTPUT}/All.public_epitope_annotations.Zscores.csv
if [ -f ${f} ] && [ ! -w ${f} ] ; then
	echo "Write-protected ${f} exists. Skipping."
else
	join --header -t, /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.join_sorted.csv \
		${OUTPUT}/All.count.Zscores.csv > ${f}
		#${OUTPUT}/All.count.Zscores.reordered.join_sorted.csv > ${f}

	chmod -w ${f}
fi







#	I think that the output of the zscore command needs to have the columns reordered.

#	They are "all of the samples ...,id,group,input"

#	Not sure if "input" is needed anymore or what "group" is.

#	I don't use this.


#awk 'BEGIN{FS=OFS=","}(NR==1){print $5,$6,$1,$2,$3,$4,$7}(NR>1){printf "%d,%d,%.2g,%.2g,%.2g,%.2g,%d\n",$5,$6,$1,$2,$3,$4,$7}' \
#  Elledge/fastq_files/merged.combined.count.Zscores.csv > Elledge/fastq_files/merged.combined.count.Zscores.reordered.csv





#	Determine actual hits by zscore threshold in both replicates

echo "Booleanizing"

while read subject ; do
	echo $subject

	all_samples=$( awk -F, -v subject=${subject} '( $1==subject && $4 != "input" ){print $2}' ${MANIFEST} | paste -sd' ' )
	echo $all_samples

	for threshold in ${THRESHOLDS}; do

		echo $threshold

		f=${OUTPUT}/${subject}.count.Zscores.${threshold}.hits.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else
			booleanize_Zscore_replicates.py --sample ${subject} \
				--threshold ${threshold} \
				--matrix ${OUTPUT}/All.count.Zscores.csv \
				--output ${f} ${all_samples}
			chmod -w ${f}
		fi

	done

	f=${OUTPUT}/${subject}.count.Zscores.minimums.csv
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		minimum_Zscore_replicates.py --sample ${subject} \
			--matrix ${OUTPUT}/All.count.Zscores.csv \
			--output ${f} ${all_samples}
		chmod -w ${f}
	fi

done < <( awk -F, '( NR>1 && $4 != "input" ){print $1}' ${MANIFEST} | sort | uniq )



f=${OUTPUT}/All.count.Zscores.minimums.csv
if [ -f ${f} ] && [ ! -w ${f} ] ; then
	echo "Write-protected ${f} exists. Skipping."
else
	#merge_all_combined_counts_files.py -o ${f} \
	merge_all_combined_counts_files.py -o ${OUTPUT}/tmp.csv \
		${OUTPUT}/*.count.Zscores.minimums.csv
	head -1 ${OUTPUT}/tmp.csv > ${f}
	tail -n +2 ${OUTPUT}/tmp.csv | sort -t, -k1,1 >> ${f}
	chmod -w ${f}
fi



if ${STOP_AFTER_ZSCORE} ; then
	exit
fi



#	**Account for public epitopes BEFORE virus score**


#	Create a list of viral species sorted by number of hits for all samples




#	process multiple Zscore thresholds

for threshold in ${THRESHOLDS}; do

	echo $threshold

	if [ -z "${SPECIES_ORDER}" ] ; then

		echo "Creating a fixed species order"

		f=${OUTPUT}/${threshold}.species_order.txt
		SPECIES_ORDER=${f}
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			tail -q -n +2 ${OUTPUT}/*.count.Zscores.${threshold}.hits.csv | sort -t, -k1,1 \
				| awk -F, '($2=="True")' > ${OUTPUT}/All.count.Zscores.${threshold}.merged_trues.csv

			sed -i '1iid,all' ${OUTPUT}/All.count.Zscores.${threshold}.merged_trues.csv

			join --header -t, ${OUTPUT}/All.count.Zscores.${threshold}.merged_trues.csv \
				/francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv > ${OUTPUT}/tmp

			awk -F, '(NR>1){print $3}' ${OUTPUT}/tmp | sort | uniq -c | sort -k1nr,1 | sed 's/^ *//' \
				| cut -d' ' -f2- > ${f}
			\rm ${OUTPUT}/tmp

			chmod -w ${f}
		fi

	fi



	#	Create virus scores in same order for all samples using the previous species order file.


	echo "Calculate virus scores"

	while read subject ; do
		echo $subject

		f=${OUTPUT}/${subject}.count.Zscores.${threshold}.hits.virus_scores.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else
			elledge_calc_scores_nofilter_forceorder.py --hits ${OUTPUT}/${subject}.count.Zscores.${threshold}.hits.csv \
				--oligo_metadata /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.csv \
				--species_order ${SPECIES_ORDER} > ${OUTPUT}/tmp 
			head -1 ${OUTPUT}/tmp > ${f}
			tail -n +2 ${OUTPUT}/tmp | sort -t, -k1,1 >> ${f}
			\rm ${OUTPUT}/tmp

			chmod -w ${f}
		fi

	done < <( awk -F, '(NR>1 && $4 != "input" ){print $1}' ${MANIFEST} | sort | uniq )



	f=${OUTPUT}/merged.${threshold}.virus_scores.csv
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		merge_results.py --int -o ${f} ${OUTPUT}/*.${threshold}.hits.virus_scores.csv
		chmod -w ${f}
	fi








	#	Threshold check


	#	~/github/ucsffrancislab/PhIP-Seq/Elledge/VirScan_viral_thresholds.csv 


	#	A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file “VirScan_viral_thresholds” contains the thresholds for each virus (Supplementary materials).


	#	GREATER THAN THRESHOLD. NOT GREATER THAN OR EQUAL TO THE THRESHOLD.


	echo "Thresholds"

	for virus_scores in ${OUTPUT}/*.count.Zscores.${threshold}.hits.virus_scores.csv ; do
		echo ${virus_scores}

		f=${virus_scores%.csv}.threshold.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else
			join --header -t, ~/github/ucsffrancislab/PhIP-Seq/Elledge/VirScan_viral_thresholds.csv ${virus_scores} \
				| awk 'BEGIN{FS=OFS=","}(NR==1){print "Species",$3}(NR>1 && $3>$2){print $1,$3}' \
				> ${f}
			chmod -w ${f}
		fi

	done 





	echo "Filter with join to public epitopes BEFORE"

	while read subject ; do
		echo $subject

		hits=${OUTPUT}/${subject}.count.Zscores.${threshold}.hits.csv

	#for hits in ${OUTPUT}/*.count.Zscores.${threshold}.hits.csv ; do
		echo ${hits}

		f=${hits%.csv}.found_public_epitopes.BEFORE_scoring.txt
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( tail -n +2 ${hits} | sort -t, -k1,1 ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '($2=="True"){print $6}' | sort | uniq > ${f}
			sed -i '1iSpecies' ${f}

			chmod -w ${f}
		fi

		f=${hits%.csv}.found_public_epitope_counts.BEFORE_scoring.txt
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( tail -n +2 ${hits} | sort -t, -k1,1 ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '($2=="True"){print $6}' | sort | uniq -c  | sort -k1nr,1 | sed -e 's/^\s*//' -e 's/ /,/' > ${f}
			sed -i '1icount,Species' ${f}

			chmod -w ${f}
		fi



		#f=${hits%.csv}.found_public_epitope_counts.BEFORE_scoring.test.csv
		f=${hits%.csv}.found_public_epitope_counts.BEFORE_scoring.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( tail -n +2 ${hits} | sort -t, -k1,1 ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '($2=="True"){print $6}' | sort | uniq -c  | sort -k1nr,1 | sed -e 's/^\s*//' -e 's/ /,/' \
				| awk 'BEGIN{FS=OFS=","}{print $2,$1}' > ${f}
			sed -i "1iSpecies,${subject}" ${f}

			chmod -w ${f}
		fi

	#done
	done < <( awk -F, '(NR>1 && $4 != "input" ){print $1}' ${MANIFEST} | sort | uniq )


	f=${OUTPUT}/merged.${threshold}.public_epitopes_BEFORE.csv
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		#merge_results.py --int -o ${f} ${OUTPUT}/*found_public_epitope_counts.BEFORE_scoring.test.csv
		merge_results.py --int -o ${f} ${OUTPUT}/*.${threshold}.hits.found_public_epitope_counts.BEFORE_scoring.csv
		chmod -w ${f}
	fi



	echo "Filter with join to public epitopes AFTER"

	#	*.7.peptides.txt files created during my version of the calc_virus_score scripts
	#	This is the list of peptides that survive the "novel" test hence the name "AFTER"
	#	The "BEFORE" list would effectively be the hits list, but I need to create something to get the counts for both BEFORE and AFTER.

	while read subject ; do
		echo $subject

		peptides=${OUTPUT}/${subject}.count.Zscores.${threshold}.hits.7.peptides.txt
	#for peptides in ${OUTPUT}/*.count.Zscores.hits.7.peptides.txt ; do
		echo ${peptides}

		f=${peptides%.7.peptides.txt}.found_public_epitopes.AFTER_scoring.txt
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( sort -t, -k1,1 ${peptides} ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '{print $2}' | sort | uniq > ${f}
			sed -i '1iSpecies' ${f}

			chmod -w ${f}
		fi

		f=${peptides%.7.peptides.txt}.found_public_epitope_counts.AFTER_scoring.txt
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( sort -t, -k1,1 ${peptides} ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '{print $2}' | sort | uniq -c  | sort -k1nr,1 | sed -e 's/^\s*//' -e 's/ /,/' > ${f}
			sed -i '1icount,Species' ${f}

			chmod -w ${f}
		fi



		#f=${peptides%.7.peptides.txt}.found_public_epitope_counts.AFTER_scoring.test.csv
		f=${peptides%.7.peptides.txt}.found_public_epitope_counts.AFTER_scoring.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else

			join -t, <( sort -t, -k1,1 ${peptides} ) \
				<( tail -n +2 /francislab/data1/refs/PhIP-Seq/public_epitope_annotations.sorted.csv ) \
				| awk -F, '{print $2}' | sort | uniq -c  | sort -k1nr,1 | sed -e 's/^\s*//' -e 's/ /,/' \
				| awk 'BEGIN{FS=OFS=","}{print $2,$1}' > ${f}
			sed -i "1iSpecies,${subject}" ${f}

			chmod -w ${f}
		fi

	#done
	done < <( awk -F, '(NR>1 && $4 != "input" ){print $1}' ${MANIFEST} | sort | uniq )



	f=${OUTPUT}/merged.public_epitopes_AFTER.csv
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		merge_results.py --int -o ${f} ${OUTPUT}/*.${threshold}.hits.7.peptides.found_public_epitope_counts.AFTER_scoring.csv
		chmod -w ${f}
	fi




	for scoring in ${OUTPUT}/*.count.Zscores.${threshold}.hits.found_public_epitopes.*_scoring.txt ; do
		echo ${scoring}

		f=${scoring%.txt}.seropositive.csv
		if [ -f ${f} ] && [ ! -w ${f} ] ; then
			echo "Write-protected ${f} exists. Skipping."
		else
			join --header -t, ${scoring} \
				${scoring%.found_public_epitopes.*_scoring.txt}.virus_scores.threshold.csv > ${f}

	#	merging gets its name from the filename, not the column name. Should I change that?
			if [[ "${scoring}" =~ "BEFORE" ]]; then
				sed -i '1s/$/_B/' ${f}
			elif [[ "${scoring}" =~ "AFTER" ]]; then
				sed -i '1s/$/_A/' ${f}
			fi

			chmod -w ${f}
		fi

	done 



	f=${OUTPUT}/merged.${threshold}.seropositive.csv
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected ${f} exists. Skipping."
	else
		#merge_results.py --int -o ${f} ${OUTPUT}/*_scoring.seropositive.csv
		merge_results.py --int -o ${f} ${OUTPUT}/*.${threshold}.hits.found_public_epitopes.*_scoring.txt
		chmod -w ${f}
	fi


done	#	thresholds






