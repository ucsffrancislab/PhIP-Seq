
#	PhIP-Seq

Phage ImmunoPrecipitation Sequencing

Collection of scripts and references on the subject from multiple sources.

Processing is predominantly based on 
[Elledge 2022 - VirScan: High-throughput Profiling of Antiviral Antibody Epitopes](https://en.bio-protocol.org/en/bpdetail?id=4464&type=0)



##	Installation

Install the custom R packages available through the paper or this repository ...
```
module load r

R -e "install.packages(c('dplyr','data.table','magrittr','R.utils','readxl','readr','VGAM','MASS','tidyr'))"
R CMD INSTALL mmR_0.1.0.tar.gz 
R CMD INSTALL virScanR_0.1.0.9000.tar.gz
```

###	Build a reference

Based on this paper, create this short sequence reference ...
```
bowtie-build vir3.fa vir3
```

However, if you wish to build references from the entire sequences, extract, trim, sort and uniq them with ...
```
zcat VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")"}(NR>1){o=substr($18,17,168);print $17","o}' | uniq | sort -t, -k1n,1 | uniq | awk -F, '{print ">"$1;print $2}' | gzip > VIR3_clean.uniq.fna.gz
```

Then build away ...
```
bowtie-build VIR3_clean.uniq.fna.gz VIR3_clean
bowtie2-build VIR3_clean.uniq.fna.gz VIR3_clean
```





##	VIR3 Post Sequencing Analysis


Scripts based on Elledge's 2022 paper's scripts or instructions.




###	Notable questions


This paper was made using a 50bp reference sequences, rather than the entire 168bp. Why?

This paper uses example single-ended reads that are 75bp then trims 25bp. Why?
Just to match the reference? 
Should we always use a reference of the exact same length of the sequences?



This paper assumes that you have SOME NUMBER of blank no-serum controls to be used as "input" in Z-score computation.
Is the specific number important?
Are they needed?
Does the R script require this specific number?
Does it need this column at all?
Is there a value that can be used in place of this?





This paper assumes that you always have 2 technical replicates. Why?
I think that you can get away with having different numbers so long as they are merged properly.

A Z-score of at least 3.5 in both technical replicates of a sample is required to call a peptide a "hit". Why 3.5?




This paper does not include all of the viruses in `VIR3_clean` in its thresholds so they are ignored. Why?

virus seropositivity thresholds? I computed something like 250:1. How were these determined? Why that way?



How do I determine "if at least one public epitope from that virus scores as a hit"?










###	Individual

####	[Align and extract counts for each FASTQ file](scripts/elledge_process_sample.bash)

* INPUT (FASTQ FILE)


This is a simple wrapper script that aligns with bowtie and extracts alignment counts from the index with `samtools idxstats`.
Using `idxstats` is fast, but you CANNOT select alignments based on quality score.


```
for f in Elledge/fastq_files/*fastq.gz ; do
  s=$( echo $(basename $f .fastq.gz)| cut -d_ -f3 )
  elledge_process_sample.bash --sample_id ${s} --index ~/github/ucsffrancislab/PhIP-Seq/Elledge/vir3 ${f} 
done
```


###	Aggregation

####	[sum all counts files for each sample](scripts/sum_counts_files.py)

* --output COMBINED_COUNTS_CSV_FILE
* List or glob of all counts files for the sample



I don't think that we've ever used had a single sample sequenced over multiple lanes,
but if we do this is how they are to be combined. 
This really just takes multiple csv files and sums up the values.



```
for f in Elledge/fastq_files/LIB0*_L001_R1.count.csv.gz ; do
  i=${f/L001/L00*}
  o=${f/L001/L001_2_3_4}
  sum_counts_files.py --output ${o} ${i}
done
```



####	[Merge all combined counts files](scripts/merge_all_combined_counts_files.py)

* --output MERGED_COUNTS_CSV_FILE
* List or glob of all combined counts files



This is tricky.
The intent is to create a count matrix that looks like ...

```
id,S148,S150,S154,S156,input
1,0,4,1,0,4
2,39,55,58,41,266
3,0,4,7,0,18
4,6,15,10,9,48
5,86,71,85,93,346
```

1 row for each virus sequence
1 column for each sample
1 final column for a summation of all of the input samples.

*NOTE: I DID NOT DEMONSTRATE HOW TO CREATE THE input count csv. The sum_counts script would likely work.*





No no-serum control sample data was provided, but the counts were, so extracting them.
```
awk 'BEGIN{FS=OFS=","}{print $1,$6}' Elledge/counts_vir3_Bio_Protocol.csv | gzip > Elledge/fastq_files/input.count.csv.gz
```

```
merge_all_combined_counts_files.py --output Elledge/fastq_files/merged.combined.count.csv \
  Elledge/fastq_files/*L001_2_3_4_R1.count.csv.gz Elledge/fastq_files/input.count.csv.gz 
```


*NOTE: The resulting table may very well need the column headers editted.*


The provided csv has a few items out or order, trailing ^M's and no EOF.
```
head -1 Elledge/counts_vir3_Bio_Protocol.csv | tr -d '\r' > Elledge/counts_vir3_Bio_Protocol.sorted.csv
tail -n +2 Elledge/counts_vir3_Bio_Protocol.csv | tr -d '\r' | sort -t, -k1n,1 >> Elledge/counts_vir3_Bio_Protocol.sorted.csv
```


Compare to example counts ...
```
diff Elledge/fastq_files/merged.combined.count.csv Elledge/counts_vir3_Bio_Protocol.sorted.csv
```
The same!



####	[Calculate Z-scores](scripts/elledge_Zscore_analysis.R)

* Merged combined counts csv file





Again, this requires no-serum controls to create a "input" column.
Not sure if a specific number of these controls are needed.
Not sure what to do if you don't have any.






```
elledge_Zscore_analysis.R Elledge/fastq_files/merged.combined.count.csv
```


Can't seem to control the column order or the number of sig-figs

```
head Elledge/fastq_files/merged.combined.count.Zscores.csv Elledge/Zscores_vir3.csv 
==> Elledge/fastq_files/merged.combined.count.Zscores.csv <==
S148,S150,S154,S156,id,group,input
-0.705073916469041,2.98925858133891,0.208099553465372,-0.777807289938226,1,5,4
-0.362792856470446,0.237458091552846,1.24173256094102,-1.06579692537454,2,217,266
-1.31561755837116,0.350848956004743,1.98886536972319,-1.36930134848176,3,19,18
-0.314599488042471,1.74107816510524,0.592221808619408,-0.265171749479516,4,49,48
1.6735612069963,0.218371255473301,1.81834269944266,0.863125456314862,5,235,346
-1.0231998371897,1.04103466642759,1.76618702852712,0.749113637034995,6,13,12
51.9983796083297,44.1481506457624,40.0726896625729,52.932337165844,7,6,5
0.163355315623413,1.87659454129098,0.468362753941683,0.800365543672742,8,32,31
-0.0340581521710463,-0.903801827362798,1.58710272421741,-0.955269414977745,9,8,7

==> Elledge/Zscores_vir3.csv <==
id,group,S148,S150,S154,S156,input
1,5,-0.71,2.99,0.21,-0.78,4
2,217,-0.36,0.24,1.24,-1.07,266
3,19,-1.32,0.35,1.99,-1.37,18
4,49,-0.31,1.74,0.59,-0.27,48
5,235,1.67,0.22,1.82,0.86,346
6,13,-1.02,1.04,1.77,0.75,12
7,6,52,44.15,40.07,52.93,5
8,32,0.16,1.88,0.47,0.8,31
9,8,-0.03,-0.9,1.59,-0.96,7
```


Geno suggests trying something like ...
```
x = data.frame(lapply(x, function(y) if(is.numeric(y)) round(y, 2) else y)) 
```
but doesn't explain why they are different.




```
awk 'BEGIN{FS=OFS=","}(NR==1){print $5,$6,$1,$2,$3,$4,$7}(NR>1){printf "%d,%d,%.2g,%.2g,%.2g,%.2g,%d\n",$5,$6,$1,$2,$3,$4,$7}' \
  Elledge/fastq_files/merged.combined.count.Zscores.csv > Elledge/fastq_files/merged.combined.count.Zscores.reordered.csv
```

```
sdiff -s Elledge/Zscores_vir3.csv Elledge/fastq_files/merged.combined.count.Zscores.reordered.csv
```

Even after reordering and rounding / trimming with printf, the numbers are not identical.
Nearly every row has at least 1 minor difference.




####	[Booleanize Zscores at threshold of 3.5 for each technical replicate](scripts/booleanize_Zscore_replicates.py)

1. List of all technical replicates
2. Z-scores file



*QUESTION : Why is a cutoff of 3.5 is used?*


The example data has only 148 and 150.

```
head Elledge/hits_combined_vir3_3.5_cutoff.csv
id,S148,S150
1,False,False
```

```
head Elledge/sample_legend_vir3.csv 
replicate_1,replicate_2,sample,protocol,IP_reagent,library,sample
S148,S154,sample_1,standard single IP,Protein AG,vir3,0.2 uL human sera
S150,S156,sample_2,standard single IP,Protein AG,vir3,0.2 uL human sera
```



```
booleanize_Zscore_replicates.py S148 S154 Elledge/fastq_files/merged.combined.count.Zscores.csv
booleanize_Zscore_replicates.py S150 S156 Elledge/fastq_files/merged.combined.count.Zscores.csv
```

####	[Merge booleanized Z-scores](scripts/merge_booleanized_replicates.py)

* --output MERGED_COMBINED_COUNTS_ZSCORE_BOOLEANIZED_REPLICATES_CSV_FILE
* List or glob of all replicates zscores combined counts files



I don't see this output file being used after this point and it may have little use.




```
merge_booleanized_replicates.py --output Elledge/fastq_files/merged.combined.count.Zscores.booleanized_replicates.csv \
  Elledge/fastq_files/merged.combined.count.Zscores.S???.csv
```

```
head -1 Elledge/hits_combined_vir3_3.5_cutoff.csv | tr -d '\r' > Elledge/hits_combined_vir3_3.5_cutoff.sorted.csv
tail -n +2 Elledge/hits_combined_vir3_3.5_cutoff.csv | tr -d '\r' | sort -t, -k1n,1 >> Elledge/hits_combined_vir3_3.5_cutoff.sorted.csv
```

Only a few differences likely due to the rounding issues 
```
sdiff -s merged.combined.count.Zscores.booleanized_replicates.csv Elledge/hits_combined_vir3_3.5_cutoff.sorted.csv
39869,True,False					      |	39869,False,False
57910,True,True						      |	57910,False,True
86690,False,True					      |	86690,False,False
89812,True,True						      |	89812,True,False
```

This output file doesn't seem to be used elsewhere or useful.






####	[Calculate virus scores](scripts/elledge_calc_scores_nofilter.py)

* hits - merged replicate zscores ( merged.combined.count.Zscores.S148.csv )
* oligo_metadata - Elledge/VIR3_clean.csv.gz ( we'll need to make one of these for our data )
* level - Species
* epitope_len - 7





Here, I think, we using the booleanize sequence zscore, a virus score is created.





```
for i in Elledge/fastq_files/merged.combined.count.Zscores.S???.csv ; do
  elledge_calc_scores_nofilter.py $i Elledge/VIR3_clean.csv.gz Species 7 > ${i%.csv}.virus_scores.csv
done
```


A few differences here.

```
sdiff -sW Elledge/virus_scores/virus_scores_S148.csv Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv
Species,SAMPLE_ID       				      |	Species,S148
Alphapapillomavirus 9,1 				      |	Alphapapillomavirus 9,2
BK polyomavirus (BKPyV),1       			      |	BK polyomavirus (BKPyV),0
Rhinovirus A,1  					      |	Rhinovirus A,2
Simian virus 12,0       				      |	Simian virus 12,1

sdiff -sW Elledge/virus_scores/virus_scores_S150.csv Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv
Species,SAMPLE_ID       				      |	Species,S150
Dengue virus,6  					      |	Dengue virus,7
```


####	Determining virus seropositivity

A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file “VirScan_viral_thresholds” contains the thresholds for each virus (Supplementary materials).
*NOTE : Public epitope annotations are available upon request.*





Not sure how to determine "if at least one public epitope from that virus scores as a hit" as I have no idea what the public epitopes are. I may have to request them as noted.




Based on the above line from the paper, I think that the following is accurate.



The virus score is the number of "novel" peptides for a given virus Species marked TRUE.
Novel determined by the longest match being less than the passed epitope_len (7)
Larger epitope_len would create higher virus scores.

This means that each virus score is dependent on the number of "novel" peptides.


This is why they use a separate threshold for each virus. 
Not sure where these treshold came from.

*NOTE : This threshold file only has thresholds for 206 viruses while the reference contains more than 440.*


*QUESTION : Why only half of the viruses in thresholds? Other irrelevant?*

*QUESTION : How were the thresholds computed? They appear to be in a ratio of 250 peptides to threshold of 1*


Join with threshold and select 

Changed "Non-A, non-B hepatitis virus" to "Non-A/non-B hepatitis virus" in both files


```
join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3)'


Alphapapillomavirus 5,2,1
Alphapapillomavirus 9,2,1
Betacoronavirus 1,4,1
Chikungunya virus,2,1
Cowpox virus,3,1.613518806
Dengue virus,8,1
Enterovirus B,5,1
Enterovirus C,2,1
Hepatitis B virus,10,1
Hepatitis C virus,8,1
Hepatitis E virus,6,1
Human adenovirus A,2,1
Human adenovirus B,2,1
Human adenovirus C,9,1
Human adenovirus F,4,1
Human coronavirus HKU1,3,1
Human herpesvirus 1,13,2.954080914
Human herpesvirus 2,11,2.633299984
Human herpesvirus 3,9,2.562639888
Human herpesvirus 4,42,2.995824612
Human herpesvirus 5,21,5.485330838
Human herpesvirus 6A,3,2.716506417
Human immunodeficiency virus 1,5,1
Human respiratory syncytial virus,4,1
Human rotavirus B219,2,1
Influenza A virus,10,1
Influenza B virus,6,1
Influenza C virus,3,1
KI polyomavirus,2,1
Lagos bat virus,2,1
Louping ill virus,2,1
Macacine herpesvirus 1,9,2.164826461
Marburg marburgvirus,3,1
Molluscum contagiosum virus,17,3.859539813
Monkeypox virus,2,1
Mumps virus,2,1
Orf virus,12,3.846913593
Papiine herpesvirus 2,11,1.933811245
Rhinovirus A,2,1
Rhinovirus B,2,1
Rotavirus A,3,1
Rotavirus B,2,1
Rubella virus,3,1
Sandfly fever Naples virus,2,1
Torque teno virus,2,1
Vaccinia virus,5,3.931113011
Yellow fever virus,3,1
```




```
join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3)'

Aichivirus A,2,1
Alphacoronavirus 1,3,1
Alphapapillomavirus 1,2,1
Alphapapillomavirus 2,2,1
Alphapapillomavirus 7,3,1
Australian bat lyssavirus,2,1
Banna virus,2,1
Betacoronavirus 1,2,1
Betapapillomavirus 1,2,1
Bunyamwera virus,2,1
Chikungunya virus,2,1
Cowpox virus,2,1.613518806
Dengue virus,7,1
Eastern equine encephalitis virus,3,1
Enterovirus B,6,1
Enterovirus C,3,1
Hepatitis B virus,12,1
Hepatitis C virus,8,1
Hepatitis E virus,5,1
Human adenovirus B,2,1
Human adenovirus C,2,1
Human adenovirus D,5,1
Human adenovirus E,3,1
Human adenovirus F,4,1
Human coronavirus HKU1,4,1
Human herpesvirus 1,29,2.954080914
Human herpesvirus 2,12,2.633299984
Human herpesvirus 3,8,2.562639888
Human herpesvirus 4,8,2.995824612
Human herpesvirus 5,18,5.485330838
Human herpesvirus 6B,3,1.344059842
Human herpesvirus 8,9,2.837298795
Human immunodeficiency virus 1,6,1
Human immunodeficiency virus 2,4,1
Human respiratory syncytial virus,5,1
Influenza A virus,17,1
Influenza B virus,7,1
Influenza C virus,2,1
Macacine herpesvirus 1,10,2.164826461
Mamastrovirus 1,3,1
Molluscum contagiosum virus,12,3.859539813
Mumps virus,2,1
Norwalk virus,4,1
Orf virus,20,3.846913593
Papiine herpesvirus 2,7,1.933811245
Parainfluenza virus 5,3,1
Primate T-lymphotropic virus 1,2,1
Pseudocowpox virus,4,1.120803941
Puumala virus,3,1
Rhinovirus A,4,1
Rhinovirus B,3,1
Rotavirus A,6,1
Rotavirus B,4,1
Rubella virus,3,1
Venezuelan equine encephalitis virus,3,1
Yellow fever virus,2,1
Zaire ebolavirus,2,1
```


Most thresholds are 1 making it have minimal effect as a filter.



