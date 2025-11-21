#!/usr/bin/env bash
#SBATCH --export=NONE   # required when using 'module'

hostname
echo "Slurm job id:${SLURM_JOBID}:"
date
start_seconds=$SECONDS



set -e	#	exit if any command fails
set -u	#	Error on usage of unset variables
set -o pipefail
if [ -n "$( declare -F module )" ] ; then
	echo "Loading required modules"
	module load CBI r
fi
#set -x	#	print expanded command before executing it
set -v


#	/francislab/data1/working/20250409-Illumina-PhIP/20250411-PhIP/out.plate5

#	/francislab/data1/working/20250409-Illumina-PhIP/20250411-PhIP/out.plate6

OUTPUT=$1


echo "Computing Zscores"

#f=${OUTPUT}/All.count.Zscores.csv
#if [ -f ${f} ] && [ ! -w ${f} ] ; then
#  echo "Write-protected ${f} exists. Skipping."
#else
#  echo "Creating ${f}"


# elledge script output is mucked up


#	#\rm ${OUTPUT}/All.count.Zscores.csv

	elledge_Zscore_analysis.R ${OUTPUT}/All.count.csv
	mv ${OUTPUT}/All.count.Zscores.csv ${OUTPUT}/tmp1.csv

	awk 'BEGIN{FS=OFS=","}{x=$(NF-2);NF=NF-3;print x,$0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv
	#awk 'BEGIN{FS=OFS=","}{x=$(NF-2);print x,$0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv

	head -1 ${OUTPUT}/tmp2.csv > ${OUTPUT}/All.count.Zscores.elledge.05.csv
	tail -n +2 ${OUTPUT}/tmp2.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.elledge.05.csv

## only 62 columns for some reason
#
##	> virScanR::vs.convert(data = counts[1:20,c(1:52, 80:90)], paramsList = convert_params1)$out
#
#
#	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp3.csv
#	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp3.csv > ${OUTPUT}/tmp4.csv
#	mv ${OUTPUT}/tmp3.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake1.csv
#
#	head -1 ${OUTPUT}/tmp4.csv > ${OUTPUT}/All.count.Zscores.jake1.csv
#	tail -n +2 ${OUTPUT}/tmp4.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake1.csv
#
#
#	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp5.csv --threshold_adjustment_factor 1.17
#	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp5.csv > ${OUTPUT}/tmp6.csv
#	mv ${OUTPUT}/tmp5.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake117.csv
#
#	head -1 ${OUTPUT}/tmp6.csv > ${OUTPUT}/All.count.Zscores.jake117.csv
#	tail -n +2 ${OUTPUT}/tmp6.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake117.csv
#
#
#	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp7.csv --threshold_adjustment_factor 1.33
#	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp7.csv > ${OUTPUT}/tmp8.csv
#	mv ${OUTPUT}/tmp7.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake133.csv
#
#	head -1 ${OUTPUT}/tmp8.csv > ${OUTPUT}/All.count.Zscores.jake133.csv
#	tail -n +2 ${OUTPUT}/tmp8.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake133.csv
#
#
#	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp9.csv --threshold_adjustment_factor 1.77
#	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp9.csv > ${OUTPUT}/tmp10.csv
#	mv ${OUTPUT}/tmp9.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake177.csv
#
#	head -1 ${OUTPUT}/tmp10.csv > ${OUTPUT}/All.count.Zscores.jake177.csv
#	tail -n +2 ${OUTPUT}/tmp10.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake177.csv




#	zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp9.csv --fixed_tile_count 300
#	awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp9.csv > ${OUTPUT}/tmp10.csv
#	mv ${OUTPUT}/tmp9.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake300.04.csv
#
#	head -1 ${OUTPUT}/tmp10.csv > ${OUTPUT}/All.count.Zscores.jake300.04.csv
#	tail -n +2 ${OUTPUT}/tmp10.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake300.04.csv
#




#  chmod -w ${f}
#  #\rm ${OUTPUT}/tmp1.csv ${OUTPUT}/tmp2.csv
#  \rm ${OUTPUT}/tmp?.csv
#
#fi



echo "Done"
date

end_seconds=$SECONDS
seconds=$((end_seconds-start_seconds))
printf "%d seconds passed\n" $seconds
echo "Elapsed: $(($seconds / 3600))hrs $((($seconds / 60) % 60))min $(($seconds % 60))sec"


