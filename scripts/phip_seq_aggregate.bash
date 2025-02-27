#!/usr/bin/env bash

set -x


manifest=$1
dir=$2



for f in ${dir}/merged.*.virus_scores.csv ; do
	echo
	echo "Annotate Virus Scores"
	echo $f
	threshold=$( basename $f .virus_scores.csv )
	threshold=${threshold#merged.}
	echo $threshold



	#head -1 ${f} | sed -e '1s/\(,[^,]*\)_./\1/g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
	#	this didn't work in plate2 and mucked things up
	#	not even sure that it was needed. id and subject are the same thing at this point!
	head -1 ${f} | sed -e '1s/\(,[^,]*\)/\1/g' -e '1s/^id/subject/' > ${dir}/tmp1.csv
	sed -e '1s/\(,[^,]*\)/\1'${i}'/g' ${f} >> ${dir}/tmp1.csv



	cat ${dir}/tmp1.csv | datamash transpose -t, > ${dir}/tmp2.csv

	head -1 ${dir}/tmp2.csv > ${dir}/tmp3.csv
	tail -n +2 ${dir}/tmp2.csv | sort -t, -k1,1 >> ${dir}/tmp3.csv

	cut -d, -f1,4,6 ${manifest} | head -1 > ${dir}/tmp4.csv
	cut -d, -f1,4,6 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp4.csv

	#join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) ${dir}/tmp3.csv > ${dir}/virus_scores.csv
	join --header -t, ${dir}/tmp4.csv ${dir}/tmp3.csv > ${dir}/virus_scores.${threshold}.csv
	cat ${dir}/virus_scores.${threshold}.csv | datamash transpose -t, > ${dir}/virus_scores.${threshold}.t.csv

done



#if [ -f $dir/All.count.Zscores.minimums.csv ] ; then
#	echo
#	echo "Creating/Annotating Zscores.minimums.filtered.t.csv"
#	echo
#	echo "- select only those tiles in All Plibs"
#	join -t, --header out.all/AllPlibs.csv $dir/All.count.Zscores.minimums.csv > ${dir}/tmp0.csv
#	head -5 ${dir}/tmp0.csv | cut -c1-100
#
#	#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	#	10007,-0.5776876416792566,-0.4664278247594301,-0.5
#
#	echo "- sort (probably not needed as is already sorted)"
#	head -1 ${dir}/tmp0.csv > ${dir}/tmp1.csv
#	tail -q -n +2 ${dir}/tmp0.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv
#	head -5 ${dir}/tmp1.csv | cut -c1-100
#
#	#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	#	10007,-0.5776876416792566,-0.4664278247594301,-0.5
#
#	echo "- transpose and sort"
#	cat ${dir}/tmp1.csv | datamash transpose -t, | head -1 > ${dir}/tmp2.csv
#	cat ${dir}/tmp1.csv | datamash transpose -t, | tail -n +2 | sort -t, -k1,1 >> ${dir}/tmp2.csv
#	head -5 ${dir}/tmp2.csv | cut -c1-100
#
#	#	id,100,1000,10006,10007,1001,10012,10015,10016,100
#	#	14061-01,-0.7610039640743376,-0.4968270548456315,-
#	#	14091-01,-0.5636348697154461,-0.4451715846162984,-
#	#	14160-01,-0.4114267078983219,-0.4717772169888741,-
#	#	14167-01,-0.5059288594504003,-0.4100912676948724,-
#
#	echo "- transpose back?"
#	cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv
#	head -5 ${dir}/tmp3.csv | cut -c1-100
#
#	#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	#	10007,-0.5776876416792566,-0.4664278247594301,-0.5
#
#	#	 I don't think this is necessary
#	#echo "- sort again?"
#	#head -2 ${dir}/tmp3.csv > ${dir}/tmp4.csv
#	#tail -n +3 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv
#
#	#	id,14061-01,14091-01,14160-01,14167-01,14178-01,14
#	#	100,-0.7610039640743376,-0.5636348697154461,-0.411
#	#	1000,-0.4968270548456315,-0.4451715846162984,-0.47
#	#	10006,-0.2754206710226481,-0.2416204063333634,-0.2
#	#	10007,-0.5776876416792566,-0.4664278247594301,-0.5
#
#	echo "- join with species"
#	#join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp4.csv > ${dir}/tmp5.csv
#	join --header -t, <( cut -d, -f1,2 /francislab/data1/refs/PhIP-Seq/VIR3_clean.virus_score.join_sorted.csv ) ${dir}/tmp3.csv > ${dir}/tmp5.csv
#	head -5 ${dir}/tmp5.csv | cut -c1-100
#
#	#	id,species,14061-01,14091-01,14160-01,14167-01,141
#	#	100,Human herpesvirus 3,-0.7610039640743376,-0.563
#	#	1000,Hepatitis B virus,-0.4968270548456315,-0.4451
#	#	10006,Human herpesvirus 8,-0.2754206710226481,-0.2
#	#	10007,Human herpesvirus 8,-0.5776876416792566,-0.4
#
#	echo "- transpose"
#	cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv
#	head -5 ${dir}/tmp6.csv | cut -c1-100
#
#	#	id,100,1000,10006,10007,1001,10012,10015,10016,100
#	#	species,Human herpesvirus 3,Hepatitis B virus,Huma
#	#	14061-01,-0.7610039640743376,-0.4968270548456315,-
#	#	14091-01,-0.5636348697154461,-0.4451715846162984,-
#	#	14160-01,-0.4114267078983219,-0.4717772169888741,-
#
#	echo "- create join file from manifest"
#	cut -d, -f1,4,6 ${manifest} | head -1 > ${dir}/tmp7.csv
#	cut -d, -f1,4,6 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
#	head -5 ${dir}/tmp7.csv | cut -c1-100
#
#	echo "- joining with manifest. EXPECTING A CERTAIN COLUMN FORMAT with at least 6 columns (1,4,6)"
#	echo -n "y,z," > ${dir}/Zscores.minimums.filtered.csv
#	head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.minimums.filtered.csv
#	#join --header -t, <( cut -d, -f1,4,6 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.filtered.csv
#	join --header -t, ${dir}/tmp7.csv <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.filtered.csv
#	head -5 ${dir}/Zscores.minimums.filtered.csv | cut -c1-100
#
#	#	y,z,id,100,1000,10006,10007,1001,10012,10015,10016
#	#	subject,type,group,Human herpesvirus 3,Hepatitis B
#	#	14061-01,glioma serum,case,-0.7610039640743376,-0.
#	#	14091-01,glioma serum,case,-0.5636348697154461,-0.
#	#	14160-01,glioma serum,case,-0.4114267078983219,-0.
#
#	echo "- transpose"
#	cat ${dir}/Zscores.minimums.filtered.csv | datamash transpose -t, > ${dir}/Zscores.minimums.filtered.t.csv
#	head -5 ${dir}/Zscores.minimums.filtered.t.csv | cut -c1-100
#
#	#	y,subject,14061-01,14091-01,14160-01,14167-01,1417
#	#	z,type,glioma serum,glioma serum,glioma serum,glio
#	#	id,group,case,case,case,case,case,case,case,case,c
#	#	100,Human herpesvirus 3,-0.7610039640743376,-0.563
#	#	1000,Hepatitis B virus,-0.4968270548456315,-0.4451
#
#fi



for f in ${dir}/merged.*.seropositive.csv ; do
	echo
	echo "Annotate Seropositivities"
	echo
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


if [ -f ${dir}/All.count.Zscores.csv ] ; then
	echo
	echo "Annotate Zscores"
	echo
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

fi



if [ -f ${dir}/All.count.Zscores.minimums.csv ] ; then

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
	cut -d, -f1,4,6 ${manifest} | head -1 > ${dir}/tmp7.csv
	cut -d, -f1,4,6 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
	#cut -d, -f1,4 ${manifest} | head -1 > ${dir}/tmp7.csv
	#cut -d, -f1,4 ${manifest} | tail -n +2 | sort -t, -k1,1 | uniq >> ${dir}/tmp7.csv
	head -5 ${dir}/tmp7.csv | cut -c1-100

	echo "- joining with the partial manifest"
	#echo -n "y," > ${dir}/Zscores.minimums.csv
	echo -n "y,z," > ${dir}/Zscores.minimums.csv
	head -1 ${dir}/tmp6.csv >> ${dir}/Zscores.minimums.csv
	#join --header -t, <( cut -d, -f1,4 ${manifest} | uniq ) <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.csv
	join --header -t, ${dir}/tmp7.csv <( tail -n +2 ${dir}/tmp6.csv ) >> ${dir}/Zscores.minimums.csv
	head -5 ${dir}/Zscores.minimums.csv | cut -c1-100

	echo "- transpose"
	cat ${dir}/Zscores.minimums.csv | datamash transpose -t, > ${dir}/Zscores.minimums.t.csv
	head -5 ${dir}/Zscores.minimums.t.csv | cut -c1-100

	#	hc out.plate1/Zscores.select.minimums.t.csv
	#	y,subject,14078-01,14118-01,14127-01,14142-01,14206-01,14223-01,14235-01,14337-01,14358-02,14415-02,
	#	z,type,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,gl
	#	id,group,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,c
	#	1,Papiine herpesvirus 2,-0.2909274899485718,-0.2078364948636183,-0.2146477992014841,-0.1362997723598
	#	10,Vaccinia virus,-0.2909274899485718,-0.2078364948636183,-0.2146477992014841,-0.1362997723598357,-0
	#	100,Human herpesvirus 3,-0.5265474573876068,-0.1892123505424088,-0.4333066813813382,-0.2725861511400
	#	1000,Hepatitis B virus,-0.2909274899485718,-0.2078364948636183,-0.2146477992014841,-0.13629977235983
	#	10000,Human herpesvirus 8,-0.5265474573876068,-0.4170756638471574,-0.2769096150469569,-0.27258615114
	#	10001,Human herpesvirus 8,-0.5265474573876068,-0.4170756638471574,-0.4333066813813382,-0.27258615114
	#	10002,Human herpesvirus 8,-0.3094106444623609,-0.3562587360229662,-0.4382612107274923,-0.26904781770

fi






#if [ -f ${dir}/Counts.csv ] ; then
	Q=40
	merge_all_combined_counts_files.py --de_nan --int -o ${dir}/tmp1.csv \
	  ${dir}/counts/*.q${Q}.count.csv.gz ${dir}/counts/input/*.q${Q}.count.csv.gz

	awk -F, '{print NF}' ${dir}/tmp1.csv | uniq
	#576

	wc -l ${dir}/tmp1.csv
	#113826 tmp1.csv

	head -1 ${dir}/tmp1.csv > ${dir}/tmp2.csv
	tail -n +2 ${dir}/tmp1.csv | sort -t, -k1,1 >> ${dir}/tmp2.csv

	cat ${dir}/tmp2.csv | datamash transpose -t, > ${dir}/tmp3.csv
	head -1 ${dir}/tmp3.csv > ${dir}/tmp4.csv
	tail -n +2 ${dir}/tmp3.csv | sort -t, -k1,1 >> ${dir}/tmp4.csv


	awk 'BEGIN{FS=OFS=","}{print $2,$1,$3,$4,$5,$6,$7,$8,$9}' ${manifest} > ${dir}/tmpman1.csv
	head -1 ${dir}/tmpman1.csv > ${dir}/tmpman2.csv
	tail -n +2 ${dir}/tmpman1.csv | sort -t, -k1,1 >> ${dir}/tmpman2.csv

	join --header -t, ${dir}/tmpman2.csv ${dir}/tmp4.csv > ${dir}/tmp5.csv



	cat ${dir}/tmp5.csv | datamash transpose -t, > ${dir}/tmp6.csv

#	head -7 ${dir}/tmp6.csv > ${dir}/tmp7.csv
#	sed -i 's/^/,,,,,,/' ${dir}/tmp7.csv
#	join --header -t, /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250205.HHV3.for_joining.csv <( tail -n +8 ${dir}/tmp6.csv ) >> ${dir}/tmp7.csv
#
#	awk -F, '{print NF}' ${dir}/tmp7.csv | uniq
#	#582
#
#	wc -l ${dir}/tmp7.csv
#	#1655 tmp7.csv
#
#	wc -l /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250205.HHV3.for_joining.csv
#	#1654 /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250205.HHV3.for_joining.csv
#
#	head -8 ${dir}/tmp7.csv > ${dir}/tmp8.csv
#	tail -n +9 ${dir}/tmp7.csv | sort -t, -k2,2 -k3,3 -k4n,4 -k5n,5 -k6n,6 >> ${dir}/tmp8.csv
#
#	wc -l ${dir}/tmp8.csv
#	#1655 tmp8.csv
#
#	awk -F, '{print NF}' ${dir}/tmp8.csv | uniq
#	#582
#
#	head -12 ${dir}/tmp8.csv | cut -c1-110
#	#,,,,,,sample,024JCM,024JCMdup,043MPL,043MPLdup,074KBP,074KBPdup,101VKC,101VKCdup,1301,1301dup,1302,1302dup,130
#	#,,,,,,subject,024JCM,024JCM,043MPL,043MPL,074KBP,074KBP,101VKC,101VKC,1301,1301,1302,1302,1303,1303,1304,1304,
#	#,,,,,,type,pemphigus serum,pemphigus serum,pemphigus serum,pemphigus serum,pemphigus serum,pemphigus serum,pem
#	#,,,,,,study,PEMS,PEMS,PEMS,PEMS,PEMS,PEMS,PEMS,PEMS,MENS,MENS,MENS,MENS,MENS,MENS,MENS,MENS,MENS,MENS,MENS,MEN
#	#,,,,,,group,Non Endemic Control,Non Endemic Control,Non Endemic Control,Non Endemic Control,Non Endemic Contro
#	#,,,,,,age,62,62,61,61,44,44,46,46,33,33,67,67,61,61,60,60,54,54,72,72,50,50,65,65,51,51,39,39,46,46,39,39,50,5
#	#,,,,,,sex,M,M,F,F,F,F,F,F,F,F,M,M,M,M,M,M,F,F,M,M,M,M,F,F,F,F,F,F,M,M,F,F,M,M,M,M,M,M,F,F,M,M,M,M,F,F,F,F,F,F,
#	#id,Species,Protein names,Version (entry),Version (sequence),start,end,14,14,13,13,14,14,13,13,13,13,13,13,13,1
#	#25844,Human herpesvirus 3,Alkaline nuclease (EC 3.1.-.-),47.0,1.0,1,56,2,0,1,0,11,7,0,27,3,2,1,3,0,0,0,3,0,4,0
#	#25848,Human herpesvirus 3,Alkaline nuclease (EC 3.1.-.-),47.0,1.0,113,168,0,0,2,0,0,0,0,0,1,4,3,8,1,0,1,1,0,18
#	#25849,Human herpesvirus 3,Alkaline nuclease (EC 3.1.-.-),47.0,1.0,141,196,0,28,93,8,0,4,0,0,13,2,43,1,25,5,19,
#	#25850,Human herpesvirus 3,Alkaline nuclease (EC 3.1.-.-),47.0,1.0,169,224,0,0,3,14,2,2,3,6,5,8,5,21,3,2,0,14,0
#
#	mv tmp8.csv HHV3.csv

#	head -7 ${dir}/tmp6.csv > ${dir}/tmp7.csv
#	sed -i 's/^/,,,,,,/' ${dir}/tmp7.csv
#	join --header -t, /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250207.for_joining.csv <( tail -n +8 ${dir}/tmp6.csv ) >> ${dir}/tmp7.csv



#	head -8 ${dir}/tmp6.csv > ${dir}/tmp7.csv
#	sed -i 's/^/,,,,,,/' ${dir}/tmp7.csv

#	#	This isn't uniq
#	join --header -t, /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250207.for_joining.csv <( tail -n +9 ${dir}/tmp6.csv ) >> ${dir}/tmp7.csv

	head -8 ${dir}/tmp6.csv > ${dir}/tmp7.csv
	sed -i 's/^/,/' ${dir}/tmp7.csv

	#	Use this but it only has id and species.
	join --header -t, /francislab/data1/refs/PhIP-Seq/VIR3_clean.id_species.uniq.csv <( tail -n +9 ${dir}/tmp6.csv ) >> ${dir}/tmp7.csv


	awk -F, '{print NF}' ${dir}/tmp7.csv | uniq
	#582

	wc -l ${dir}/tmp7.csv
	#126324 tmp7.csv

#	don't use this file. It isn't unique with all of the additional fields.
#	wc -l /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250207.for_joining.csv
#	#128258 /francislab/data1/refs/PhIP-Seq/VIR3_clean.20250207.for_joining.csv

	wc -l /francislab/data1/refs/PhIP-Seq/VIR3_clean.id_species.uniq.csv

	#	lose almost 2000


	head -9 ${dir}/tmp7.csv > ${dir}/Counts.csv
	#tail -n +10 ${dir}/tmp7.csv | sort -t, -k2,2 -k3,3 -k4n,4 -k5n,5 -k6n,6 >> ${dir}/Counts.csv
	tail -n +10 ${dir}/tmp7.csv | sort -t, -k2,2  >> ${dir}/Counts.csv

	#head -8 ${dir}/tmp7.csv > ${dir}/Counts.csv
	#tail -n +9 ${dir}/tmp7.csv | sort -t, -k2,2 -k3,3 -k4n,4 -k5n,5 -k6n,6 >> ${dir}/Counts.csv

#fi


phip_seq_normalize_counts.py -i ${dir}/Counts.csv
	
	
#	Prep for usage in Multi_Plate_Case_Control_Peptide_Regression.R like Zscores.csv

if [ -f ${dir}/Counts.normalized.subtracted.csv ] ; then

#	is the
#	y,x	#	probably not
#	subject,type	#	yes
#	needed?

#hc out.plate1/Zscores.csv
#y,x,id,1,10,100,1000,10000,10001,10002,10003,10004,10005,10006,10007,10008,10009,1001,10010,10011,10
#subject,type,species,Papiine herpesvirus 2,Vaccinia virus,Human herpesvirus 3,Hepatitis B virus,Huma
#14078-01,glioma serum,14078-01,-0.1735488911468417,-0.1735488911468417,-0.20851375718739676,-0.17354
#14078-01,glioma serum,14078-01dup,-0.2909274899485718,-0.2909274899485718,-0.5265474573876068,-0.290
#14118-01,glioma serum,14118-01,-0.20783649486361835,-0.20783649486361835,-0.18921235054240884,-0.207
#14118-01,glioma serum,14118-01dup,-0.1909523835530943,-0.1909523835530943,0.4679949246045931,-0.1909
#14127-01,glioma serum,14127-01,-0.2146477992014841,-0.2146477992014841,-0.4333066813813382,-0.214647
#14127-01,glioma serum,14127-01dup,-0.11776590579451997,-0.11776590579451997,0.2664847297296404,-0.11
#14142-01,glioma serum,14142-01,-0.13629977235983573,-0.13629977235983573,-0.27258615114003143,-0.136
#14142-01,glioma serum,14142-01dup,,,-0.24344136937867783,,-0.24344136937867783,-0.24344136937867783,
#
#[gwendt@c4-dev3 /francislab/data1/working/20250128-Illumina-PhIP/20250128c-PhIP]$ hc out.plate1/Counts.normalized.subtracted.csv
#id,species,protein,start,14078-01,14078-01dup,14118-01,14118-01dup,14127-01,14127-01dup,14142-01,141
#,,,,14078-01,14078-01,14118-01,14118-01,14127-01,14127-01,14142-01,14142-01,14206-01,14206-01,14223-
#,,,,/francislab/data1/working/20241204-Illumina-PhIP/20241204b-bowtie2/out/S83.VIR3_clean.1-84.bam,/
#,,,,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,glioma serum,gliom
#,,,,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,IPS,
#,,,,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,case,c
#,,,,60,60,61,61,63,63,78,78,49,49,58,58,72,72,58,58,61,61,62,62,71,71,53,53,54,54,75,75,61,61,53,53,
#,,,,M,M,M,M,M,M,M,M,M,M,M,M,M,M,F,F,F,F,M,M,M,M,M,M,M,M,F,F,F,F,F,F,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,
#,,,,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#115701,Absidia glauca (Pin mould),Actin-2,1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,

	head -2 ${dir}/Counts.normalized.subtracted.csv | tail -1 > ${dir}/tmp1.csv
	head -4 ${dir}/Counts.normalized.subtracted.csv | tail -1 >> ${dir}/tmp1.csv
	head -1 ${dir}/Counts.normalized.subtracted.csv | tail -1 >> ${dir}/tmp1.csv
	tail -n +10 ${dir}/Counts.normalized.subtracted.csv | sort -t, -k1,1 >> ${dir}/tmp1.csv

	#cat ${dir}/tmp1.csv | cut -d, -f1,2,5- | datamash transpose -t, > ${dir}/Counts.normalized.subtracted.trim.csv
	cat ${dir}/tmp1.csv | datamash transpose -t, > ${dir}/Counts.normalized.subtracted.trim.csv

#	y,x	#	probably not
#	subject,type	#	yes

	sed -i -e '1s/^,/y,x/' -e '2s/^,/subject,type/' ${dir}/Counts.normalized.subtracted.trim.csv

#	create minimums? drop the sample id. 0 out negatives

fi






echo "Cleanup"
\rm ${dir}/tmp*.csv

echo "Done"



