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


#

  elledge_Zscore_analysis.R ${OUTPUT}/All.count.csv
  mv ${OUTPUT}/All.count.Zscores.csv ${OUTPUT}/tmp1.csv
  awk 'BEGIN{FS=OFS=","}{x=$(NF-2);NF=NF-3;print x,$0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv

  head -1 ${OUTPUT}/tmp2.csv > ${OUTPUT}/All.count.Zscores.elledge.csv
  tail -n +2 ${OUTPUT}/tmp2.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.elledge.csv



  zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp3.csv #--count_threshold 200
  awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp3.csv > ${OUTPUT}/tmp4.csv
  mv ${OUTPUT}/tmp3.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake.csv

  head -1 ${OUTPUT}/tmp4.csv > ${OUTPUT}/All.count.Zscores.jake.csv
  tail -n +2 ${OUTPUT}/tmp4.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake.csv


#  chmod -w ${f}
#  #\rm ${OUTPUT}/tmp1.csv ${OUTPUT}/tmp2.csv
#  \rm ${OUTPUT}/tmp2.csv
#
#fi



echo "Done"
date

end_seconds=$SECONDS
seconds=$((end_seconds-start_seconds))
printf "%d seconds passed\n" $seconds
echo "Elapsed: $(($seconds / 3600))hrs $((($seconds / 60) % 60))min $(($seconds % 60))sec"


