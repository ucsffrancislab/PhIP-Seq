#!/usr/bin/env bash

#	/francislab/data1/working/20250409-Illumina-PhIP/20251113-PhIP/out.plate5

#	/francislab/data1/working/20250409-Illumina-PhIP/20251113-PhIP/out.plate6

OUTPUT=$1


echo "Computing Zscores"

#f=${OUTPUT}/All.count.Zscores.csv
#if [ -f ${f} ] && [ ! -w ${f} ] ; then
#  echo "Write-protected ${f} exists. Skipping."
#else
#  echo "Creating ${f}"

  elledge_Zscore_analysis.R ${OUTPUT}/All.count.csv
  mv ${OUTPUT}/All.count.Zscores.csv ${OUTPUT}/tmp1.csv
  awk 'BEGIN{FS=OFS=","}{x=$(NF-2);NF=NF-3;print x,$0}' ${OUTPUT}/tmp1.csv > ${OUTPUT}/tmp2.csv

  head -1 ${OUTPUT}/tmp2.csv > ${OUTPUT}/All.count.Zscores.elledge.csv
  tail -n +2 ${OUTPUT}/tmp2.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.elledge.csv



  zscoring.py --input ${OUTPUT}/All.count.csv --output ${OUTPUT}/tmp3.csv #--count_threshold 200
  awk 'BEGIN{FS=OFS=","}{NF=NF-2;print $0}' ${OUTPUT}/tmp3.csv > ${OUTPUT}/tmp4.csv
  mv ${OUTPUT}/tmp3.csv ${OUTPUT}/All.count.Zscores.with_input_and_bin.jake.csv

  head -1 ${OUTPUT}/tmp4.csv > ${OUTPUT}/All.count.Zscores.csv
  tail -n +2 ${OUTPUT}/tmp4.csv | sort -t, -k1,1 >> ${OUTPUT}/All.count.Zscores.jake.csv


#  chmod -w ${f}
#  #\rm ${OUTPUT}/tmp1.csv ${OUTPUT}/tmp2.csv
#  \rm ${OUTPUT}/tmp2.csv
#
#fi


