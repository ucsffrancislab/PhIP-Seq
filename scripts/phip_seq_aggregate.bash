#!/usr/bin/env bash

set -x


manifest=$1
dir=$2


echo "Annotating Zscores.minimums.filtered.t.csv"

join -t, --header All4Plibs.csv $dir/All.count.Zscores.minimums.csv > ${dir}/tmp0.csv

head -1 ${dir}/tmp0.csv > ${dir}/tmp1.csv
tail -q -n +2 ${dir}/tmp0.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv

cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv

cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv

head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv

join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp4.csv > ${dir}/tmp5.csv

cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv

echo -n "y,z," > ${dir}/Zscores.minimums.filtered.csv
head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.minimums.filtered.csv
join --header -t, <( cut -d, -f1,4,6 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.filtered.csv

cat ${dir}/Zscores.minimums.filtered.csv | datamash transpose -t, > ${dir}/Zscores.minimums.filtered.t.csv


echo "Annotate Seropositivity"

head -1 ${dir}/merged.seropositive.csv | sed -e '1s/\(,[^,]*\)_./\1/g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
sed -e '1s/\(,[^,]*\)/\1'${i}'/g' ${dir}/merged.seropositive.csv >> ${dir}/tmp1.csv

cat ${dir}/tmp1.csv | datamash transpose -t, > ${dir}/tmp2.csv

head -1 ${dir}/tmp2.csv > ${dir}/tmp3.csv
tail -n +2 ${dir}/tmp2.csv | sort -t, -k1,1 >> ${dir}/tmp3.csv

join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) ${dir}/tmp3.csv > ${dir}/seropositive.csv
cat ${dir}/seropositive.csv | datamash transpose -t, > ${dir}/seropositive.t.csv



echo "Annotate Zscores"

head -1 ${dir}/All.count.Zscores.csv | sed -e '1s/dup//g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
head -1 ${dir}/All.count.Zscores.csv >> ${dir}/tmp1.csv
tail -q -n +2 ${dir}/All.count.Zscores.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv
cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv

#	probably useless remnant
head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv
##	tmp3 == tmp4!

echo -n "x," > ${dir}/tmp5.csv
head -1 ${dir}/tmp4.csv >> ${dir}/tmp5.csv
join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) <( tail -n +2 ${dir}/tmp4.csv ) >> ${dir}/tmp5.csv

cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv

echo -n "y," > ${dir}/Zscores.csv
head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.csv
join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.csv

cat ${dir}/Zscores.csv | datamash transpose -t, > ${dir}/Zscores.t.csv

\rm ${dir}/tmp?.csv



