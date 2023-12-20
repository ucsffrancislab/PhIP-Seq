
#	PhIP-Seq

phage immunoprecipitation sequencing

Collection of scripts and references on the subject from multiple sources.




##	Methods


###	Pre

[Script based on Larman paper's instructions](scripts/larman_prep_sequences.bash)

Script takes 3 required positional arguments

1. TILESIZE (56)
2. OVERLAP (28)
3. INPUT (PATH TO FASTA FILE CONTAINING ALL AMINO ACID SEQUENCES)


###	Post


Scripts based on Elledge paper's scripts or instructions.


####	Individual

#####	[Align and extract counts for each FASTQ file](scripts/elledge_process_sample.bash)

1. INPUT (FASTQ FILE)

```
elledge_process_sample.bash Elledge/fastq_files/LIB044174_GEN00168483_S148_L001_R1.fastq.gz 
elledge_process_sample.bash Elledge/fastq_files/LIB044174_GEN00168483_S148_L002_R1.fastq.gz 
elledge_process_sample.bash Elledge/fastq_files/LIB044174_GEN00168483_S148_L003_R1.fastq.gz 
elledge_process_sample.bash Elledge/fastq_files/LIB044174_GEN00168483_S148_L004_R1.fastq.gz 
```

```
for f in Elledge/fastq_files/*fastq.gz ; do
elledge_process_sample.bash ${f}
done
```


####	Aggregation

#####	[sum all counts files for each sample](scripts/sum_counts_files.py)

* --output COMBINED_COUNTS_CSV_FILE

1. List or glob of all counts files for the sample


```
sum_counts_files.py --output Elledge/fastq_files/LIB044174_GEN00168483_S148_L001_2_3_4_R1.count.combined.csv.gz Elledge/fastq_files/LIB044174_GEN00168483_S148_L00*_R1.count.csv.gz`
```

```
for f in Elledge/fastq_files/LIB0*_L001_R1.count.csv.gz ; do
i=${f/L001/L00*}
o=${f/L001/L001_2_3_4}
sum_counts_files.py --output ${o} ${i}
done
```



#####	[Merge all combined counts files](scripts/merge_all_combined_counts_files.py)

* --output MERGED_COUNTS_CSV_FILE

1. List or glob of all combined counts files


```
merge_all_combined_counts_files.py Elledge/count.combined_files/LIB044174_GEN001684*.count.combined.csv.gz
```



No no-serum control sample data was provided, but the counts were so extracting them.
```
awk 'BEGIN{FS=OFS=","}{print $1,$6}' Elledge/counts_vir3_Bio_Protocol.csv | gzip > Elledge/fastq_files/input.count.csv.gz
```

Trimming sample names to just the S###.
```
for f in Elledge/fastq_files/*L001_2_3_4_R1.count.csv.gz ; do
l=${f/_L001_2_3_4_R1/}
l=${l/LIB*_GEN*_S/S}
ln -s $(basename $f) $l
done
```

```
merge_all_combined_counts_files.py --output Elledge/fastq_files/merged.combined.count.csv Elledge/fastq_files/S???.count.csv.gz Elledge/fastq_files/input.count.csv.gz 
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



#####	[Calculate Z-scores](scripts/elledge_Zscore_analysis.R)

1. Merged combined counts csv file


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




#####	[Booleanize Zscores at threshold of 3.5 for each technical replicate](scripts/booleanize_Zscore_replicates.py)

1. List of all technical replicate
2. Z-scores file



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

#####	[Merge booleanized Z-scores](scripts/merge_booleanized_replicates.py)

```
merge_booleanized_replicates.py --output merged.combined.count.Zscores.booleanized_replicates.csv Elledge/fastq_files/merged.combined.count.Zscores.S???.csv
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


#####	[Calculate virus scores](scripts/elledge_calc_scores_nofilter.py)

1. Booleanized Z-scores csv file


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






#####	Determining virus seropositivity

A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file “VirScan_viral_thresholds” contains the thresholds for each virus (Supplementary materials).
Note: Public epitope annotations are available upon request.

Based on the above, I think that the following is accurate.

```
awk -F, '((NR==1)||($2>3.5))' Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv
awk -F, '((NR==1)||($2>3.5))' Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv
```




##	References

###	Elledge
 
VirScan: High-throughput Profiling of Antiviral Antibody Epitopes

https://pubmed.ncbi.nlm.nih.gov/35937932/

https://doi.org/10.21769/BioProtoc.4464

https://os.bio-protocol.org/attached/file/20220629/Supplementary%20materials.docx

1) Videos: https://www.dropbox.com/sh/enlvqmsl1971bv3/AABfNDV21uO_bboLgWy3XBi4a?dl=0
2) Analysis scripts and input files: https://www.dropbox.com/sh/qvo1t75sgsq7fi8/AAAY-LQEQDrxV6wWF6OJDHPWa?dl=0
3) Example data: https://www.dropbox.com/sh/7hnqnx4yiskabgo/AAAeRE70EScRIzusqGX2HndUa?dl=0
4) Liquid handling robot protocol: https://www.dropbox.com/sh/guguiqgrmviqk9u/AABPoVam-Xsgg82su-IibDS_a?dl=0


Shrock EL, Shrock CL, Elledge SJ. VirScan: High-throughput Profiling of Antiviral Antibody Epitopes. Bio Protoc. 2022 Jul 5;12(13):e4464. doi: 10.21769/BioProtoc.4464. PMID: 35937932; PMCID: PMC9303818.

This repo includes these analysis script, input files and example data.


###	Larman

PhIP-Seq characterization of serum antibodies using oligonucleotide-encoded peptidomes

https://www.nature.com/articles/s41596-018-0025-6

https://doi.org/10.1038/s41596-018-0025-6

Nat Protoc. 2018 September ; 13(9): 1958–1978. doi:10.1038/s41596-018-0025-6

https://github.com/lasersonlab/phip-stat

https://github.com/lasersonlab/pepsyn

The Larman paper did not include any specific code or data that isn't already easily accessible.

Mohan, D., Wansley, D.L., Sie, B.M. et al. PhIP-Seq characterization of serum antibodies using oligonucleotide-encoded peptidomes. Nat Protoc 13, 1958–1978 (2018). https://doi.org/10.1038/s41596-018-0025-6


