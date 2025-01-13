#!/usr/bin/env bash

set -x


manifest=$1
dir=$2







echo "Creating/Annotating Zscores.minimums.filtered.t.csv"

echo "- select only those tiles in All Plibs"
join -t, --header out.all/AllPlibs.csv $dir/All.count.Zscores.minimums.csv > ${dir}/tmp0.csv
head -5 ${dir}/tmp0.csv | cut -c1-100

#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	10007,-0.5776876416792566,-0.4664278247594301,-0.5

echo "- sort (probably not needed as is already sorted)"
head -1 ${dir}/tmp0.csv > ${dir}/tmp1.csv
tail -q -n +2 ${dir}/tmp0.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv
head -5 ${dir}/tmp1.csv | cut -c1-100

#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	10007,-0.5776876416792566,-0.4664278247594301,-0.5

echo "- transpose and sort"
cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv
head -5 ${dir}/tmp2.csv | cut -c1-100

#	id,100,1000,10006,10007,1001,10012,10015,10016,100
#	14061-01,-0.7610039640743376,-0.4968270548456315,-
#	14091-01,-0.5636348697154461,-0.4451715846162984,-
#	14160-01,-0.4114267078983219,-0.4717772169888741,-
#	14167-01,-0.5059288594504003,-0.4100912676948724,-

echo "- transpose back?"
cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv
head -5 ${dir}/tmp3.csv | cut -c1-100

#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	10007,-0.5776876416792566,-0.4664278247594301,-0.5

#	 I don't think this is necessary
#echo "- sort again?"
#head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
#tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv

#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	10007,-0.5776876416792566,-0.4664278247594301,-0.5

echo "- join with species"
#join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp4.csv > ${dir}/tmp5.csv
join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp3.csv > ${dir}/tmp5.csv
head -5 ${dir}/tmp5.csv | cut -c1-100

#	id,species,14061-01,14091-01,14160-01,14167-01,141
#	100,Human herpesvirus 3,-0.7610039640743376,-0.563
#	1000,Hepatitis B virus,-0.4968270548456315,-0.4451
#	10006,Human herpesvirus 8,-0.2754206710226481,-0.2
#	10007,Human herpesvirus 8,-0.5776876416792566,-0.4

echo "- transpose"
cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv
head -5 ${dir}/tmp6.csv | cut -c1-100

#	id,100,1000,10006,10007,1001,10012,10015,10016,100
#	species,Human herpesvirus 3,Hepatitis B virus,Huma
#	14061-01,-0.7610039640743376,-0.4968270548456315,-
#	14091-01,-0.5636348697154461,-0.4451715846162984,-
#	14160-01,-0.4114267078983219,-0.4717772169888741,-

echo "- create join file from manifest"
cut -d, -f1,4,6 ${manifest} | head -1 > ${dir}/tmp7.csv
cut -d, -f1,4,6 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
head -5 ${dir}/tmp7.csv | cut -c1-100

echo "- joining with manifest. EXPECTING A CERTAIN COLUMN FORMAT with at least 6 columns (1,4,6)"
echo -n "y,z," > ${dir}/Zscores.minimums.filtered.csv
head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.minimums.filtered.csv
#join --header -t, <( cut -d, -f1,4,6 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.filtered.csv
join --header -t, ${dir}/tmp7.csv <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.filtered.csv
head -5 ${dir}/Zscores.minimums.filtered.csv | cut -c1-100

#	y,z,id,100,1000,10006,10007,1001,10012,10015,10016
#	subject,type,group,Human herpesvirus 3,Hepatitis B
#	14061-01,glioma serum,case,-0.7610039640743376,-0.
#	14091-01,glioma serum,case,-0.5636348697154461,-0.
#	14160-01,glioma serum,case,-0.4114267078983219,-0.

echo "- transpose"
cat ${dir}/Zscores.minimums.filtered.csv | datamash transpose -t, > ${dir}/Zscores.minimums.filtered.t.csv
head -5 ${dir}/Zscores.minimums.filtered.t.csv | cut -c1-100

#	y,subject,14061-01,14091-01,14160-01,14167-01,1417
#	z,type,glioma serum,glioma serum,glioma serum,glio
#	id,group,case,case,case,case,case,case,case,case,c
#	100,Human herpesvirus 3,-0.7610039640743376,-0.563
#	1000,Hepatitis B virus,-0.4968270548456315,-0.4451








echo
echo "Annotate Seropositivities"

for f in ${dir}/merged.*.seropositive.csv ; do
	echo $f
	threshold=$( basename $f .seropositive.csv )
	threshold=${threshold#merged.}
	echo $threshold

	head -1 ${f} | sed -e '1s/\(,[^,]*\)_./\1/g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
	sed -e '1s/\(,[^,]*\)/\1'${i}'/g' ${f} >> ${dir}/tmp1.csv

	cat ${dir}/tmp1.csv | datamash transpose -t, > ${dir}/tmp2.csv

	head -1 ${dir}/tmp2.csv > ${dir}/tmp3.csv
	tail -n +2 ${dir}/tmp2.csv | sort -t, -k1,1 >> ${dir}/tmp3.csv

	cut -d, -f1,4 ${manifest} | head -1 > ${dir}/tmp4.csv
	cut -d, -f1,4 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp4.csv

	#join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) ${dir}/tmp3.csv > ${dir}/seropositive.csv
	join --header -t, ${dir}/tmp4.csv ${dir}/tmp3.csv > ${dir}/seropositive.${threshold}.csv
	cat ${dir}/seropositive.${threshold}.csv | datamash transpose -t, > ${dir}/seropositive.${threshold}.t.csv

done






echo
echo "Annotate Zscores"

echo "- sorting"
head -1 ${dir}/All.count.Zscores.csv | sed -e '1s/dup//g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
head -1 ${dir}/All.count.Zscores.csv >> ${dir}/tmp1.csv
tail -q -n +2 ${dir}/All.count.Zscores.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv
head -5 ${dir}/tmp1.csv | cut -c1-100

echo "- transpose and sort"
cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv
head -5 ${dir}/tmp2.csv | cut -c1-100

echo "- transpose"
cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv
head -5 ${dir}/tmp3.csv | cut -c1-100

##	probably useless remnant
#head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
#tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv
###	tmp3 == tmp4!

echo "- join with virus species"
echo -n "x," > ${dir}/tmp5.csv
#head -1 ${dir}/tmp4.csv >> ${dir}/tmp5.csv
#join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) <( tail -n +2 ${dir}/tmp4.csv ) >> ${dir}/tmp5.csv
head -1 ${dir}/tmp3.csv >> ${dir}/tmp5.csv
join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) <( tail -n +2 ${dir}/tmp3.csv ) >> ${dir}/tmp5.csv
head -5 ${dir}/tmp5.csv | cut -c1-100

echo "- transpose"
cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv
head -5 ${dir}/tmp6.csv | cut -c1-100

echo "- create partial manifest join file"
cut -d, -f1,4 ${manifest} | head -1 > ${dir}/tmp7.csv
cut -d, -f1,4 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
head -5 ${dir}/tmp7.csv | cut -c1-100

echo "- join with partial manifest"
echo -n "y," > ${dir}/Zscores.csv
head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.csv
#join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.csv
join --header -t, ${dir}/tmp7.csv <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.csv
head -5 ${dir}/Zscores.csv | cut -c1-100

echo "- transpose"
cat ${dir}/Zscores.csv | datamash transpose -t, > ${dir}/Zscores.t.csv
head -5 ${dir}/Zscores.t.csv | cut -c1-100










echo
echo "Annotate Zscore minimums"

echo "- Sort All.count.Zscores.minimums.csv by tile, but probably already sorted"
#head -1 ${dir}/All.count.Zscores.minimums.csv | sed -e '1s/dup//g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
head -1 ${dir}/All.count.Zscores.minimums.csv > ${dir}/tmp1.csv
tail -q -n +2 ${dir}/All.count.Zscores.minimums.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv
head -5 ${dir}/tmp1.csv | cut -c1-100

echo "- transpose and sort by subject, but also probably already sorted"
cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv
head -5 ${dir}/tmp2.csv | cut -c1-100

echo "- tranpose back to original"
cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv
head -5 ${dir}/tmp3.csv | cut -c1-100

##	probably useless remnant
#head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
#tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv
###	tmp3 == tmp4!

#	tmp1 usually the same as tmp3

echo "- join with virus species"

#echo -n "x," > ${dir}/tmp5.csv
#head -1 ${dir}/tmp4.csv >> ${dir}/tmp5.csv
#join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) <( tail -n +2 ${dir}/tmp4.csv ) >> ${dir}/tmp5.csv
#head -1 ${dir}/tmp3.csv >> ${dir}/tmp5.csv
#join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) <( tail -n +2 ${dir}/tmp3.csv ) >> ${dir}/tmp5.csv

join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp3.csv > ${dir}/tmp5.csv
head -5 ${dir}/tmp5.csv | cut -c1-100

echo "- transpose"
cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv
head -5 ${dir}/tmp5.csv | cut -c1-100

echo "- create join file from manifest"
cut -d, -f1,4 ${manifest} | head -1 > ${dir}/tmp7.csv
cut -d, -f1,4 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
head -5 ${dir}/tmp7.csv | cut -c1-100

echo "- joining with the partial manifest"
echo -n "y," > ${dir}/Zscores.minimums.csv
head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.minimums.csv
#join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.csv
join --header -t, ${dir}/tmp7.csv <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.csv
head -5 ${dir}/Zscores.minimums.csv | cut -c1-100

echo "- transpose"
cat ${dir}/Zscores.minimums.csv | datamash transpose -t, > ${dir}/Zscores.minimums.t.csv
head -5 ${dir}/Zscores.minimums.t.csv | cut -c1-100



echo "Cleanup"
\rm ${dir}/tmp?.csv

echo "Done"



