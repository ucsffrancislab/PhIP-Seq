
#	PhIP-Seq

Phage ImmunoPrecipitation Sequencing

Collection of scripts and references on the subject from multiple sources.




##	Methods

```
module load r

R CMD INSTALL mmR_0.1.0.tar.gz 
R -e "install.packages('VGAM')"
R CMD INSTALL virScanR_0.1.0.9000.tar.gz
```


###	Prepare sequences to be synthesized

####	[Script based on Larman paper's instructions](scripts/larman_prep_sequences.bash)

Script takes 3 required positional arguments

1. TILESIZE (56)
2. OVERLAP (28)
3. INPUT (PATH TO FASTA FILE CONTAINING ALL AMINO ACID SEQUENCES)


Send the `oligos-56-28.fasta` file to a DNA synthesis company for manufacture of the oligonucleotide library.


#### Create reference

Larman creates a reference from entire oligo. Elledge trims it to just the first 50bp.

*QUESTION : Which should we do? Test them both?*

*QUESTION : Why bowtie and not bowtie2 or bwa?*

Elledge's post processing uses numbered reference reads whose metadata is in a csv.


The oligo fasta file uses reads with just a simple number as the name.
```
>1
atgcgcagcttgctgtttgtggtcggtgcttgggtcgctgctctcgtcac
>2
actacaaccaccgctgccgcagggaacacatctgcaacagcttctccagg
>3
attaccgctgccgctcctccaggtcattcaacaccttggcctgcactccc
```

This is the "id" which is carried through the whole process.

These 5 columns are referenced in the virus score calculation at the end.
Not sure if all of them are used or needed.
The default output only shows the Species, as that is a given parameter.

```
zcat Elledge/VIR3_clean.csv.gz | awk 'BEGIN{FS=OFS=","}{print $17,$12,$9,$4,$21}' | head
id,Species,Organism,Entry,peptide
1,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A126,MRSLLFVVGAWVAALVTNLTPDAALASGTTTTAAAGNTSATASPGDNATSIDAGST
2,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A126,TTTTAAAGNTSATASPGDNATSIDAGSTITAAAPPGHSTPWPALPTDLALPLVIGG
3,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A126,ITAAAPPGHSTPWPALPTDLALPLVIGGLCALTLAAMGAGALLHRCCRRCARRRQN
4,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A126,LCALTLAAMGAGALLHRCCRRCARRRQNVSSVSA
5,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A130,RDRGPSRSRVRYTRLAASEA
6,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A132,MGFGAAAALLALAVALARVPAGGGAYVPVDRALTRVSPNRFRGSSLPPPEQKTDPP
7,Papiine herpesvirus 2,Cercopithecine herpesvirus 16 (CeHV-16) (Herpesvirus papio 2),A0A132,VDRALTRVSPNRFRGSSLPPPEQKTDPPDVRRVYH
8,Vaccinia virus,Vaccinia virus,A0ERK8,MHYPKYYINITKINPHLANQFRAWKKRIAGRDYITNLSKDTGIQQSKLTETIRNCQ
9,Vaccinia virus,Vaccinia virus,A0ERK8,AGRDYITNLSKDTGIQQSKLTETIRNCQKNRNIYGLYIHYNLVINWITDVIINQY
```






#####	Our



######	Create Our Metadata file ( Elledge method )


Since we're using Elledge's technique, we'll need to do this.

Move the read and read name into a metadata file and replace it with just a number

Build a metadata file from the fasta.


```
~/github/ucsffrancislab/genomics/refs/phipSeq-20221116/oligos-ref-56-28.fasta 
>NP_000468.1_albumin_preproprotein|0-56
ATGAAATGGGTTACATTTATTAGTTTGCTTTTCCTGTTCTCTTCAGCGTATTCTCGTGGGGTTTTCCGTCGTGATGCACATAAATCAGAGGTAGCCCATCGCTTCAAGGATTTAGGGGAAGAAAACTTTAAAGCGCTGGTTCTGATTGCGTTTGCCCAATATTTGCAG
```


```
pepsyn translate ~/github/ucsffrancislab/genomics/refs/phipSeq-20221116/oligos-ref-56-28.fasta ~/github/ucsffrancislab/genomics/refs/phipSeq-20221116/oligos-ref-56-28.faa

echo "id,accession,Species,Organism,Entry,peptide" > our_vir3.csv
id=0
while read name aa ; do
id=$[id+1]
name=${name#>}
accession=$( echo $name | cut -d. -f1 )
region=$( echo $name | cut -d\| -f2- )
desc=$( echo $name | cut -d\| -f1 | cut -d_ -f3- )
echo ${id},${accession},${desc},${desc},${region},${aa}
done < <( cat ~/github/ucsffrancislab/genomics/refs/phipSeq-20221116/oligos-ref-56-28.faa | paste - - ) >> our_vir3.csv

head our_vir3.csv

id,accession,Species,Organism,Entry,peptide
1,NP_000468,albumin_preproprotein,albumin_preproprotein,0-56,MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQ
2,NP_000468,albumin_preproprotein,albumin_preproprotein,28-84,SEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAE
3,NP_000468,albumin_preproprotein,albumin_preproprotein,56-112,QCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMA
4,NP_000468,albumin_preproprotein,albumin_preproprotein,84-140,NCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLV
5,NP_000468,albumin_preproprotein,albumin_preproprotein,112-168,DCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIAR
6,NP_000468,albumin_preproprotein,albumin_preproprotein,140-196,RPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAA
7,NP_000468,albumin_preproprotein,albumin_preproprotein,168-224,RHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKC
8,NP_000468,albumin_preproprotein,albumin_preproprotein,196-252,DKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEF
9,NP_000468,albumin_preproprotein,albumin_preproprotein,224-280,ASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADD

sed -i -e 's/albumin_preproprotein/Human/' -e 's/envelope_glycoprotein_E_Human_alphaherpesvirus_3/VZV/' -e 's/Envelope_surface_glycoprotein_gp120_Human_immunodeficiency_virus_1/HIV/' our_vir3.csv

head our_vir3.csv
id,accession,Species,Organism,Entry,peptide
1,NP_000468,Human,albumin_preproprotein,0-56,MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQ
2,NP_000468,Human,albumin_preproprotein,28-84,SEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAE
3,NP_000468,Human,albumin_preproprotein,56-112,QCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMA
4,NP_000468,Human,albumin_preproprotein,84-140,NCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLV
5,NP_000468,Human,albumin_preproprotein,112-168,DCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIAR
6,NP_000468,Human,albumin_preproprotein,140-196,RPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAA
7,NP_000468,Human,albumin_preproprotein,168-224,RHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKC
8,NP_000468,Human,albumin_preproprotein,196-252,DKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEF
9,NP_000468,Human,albumin_preproprotein,224-280,ASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADD
```


```
cat ~/github/ucsffrancislab/genomics/refs/phipSeq-20221116/oligos-ref-56-28.fasta | paste - - | awk '{print ">"NR;print $2}' > our_vir3.fna

head our_vir3.fna

>1
ATGAAATGGGTTACATTTATTAGTTTGCTTTTCCTGTTCTCTTCAGCGTATTCTCGTGGGGTTTTCCGTCGTGATGCACATAAATCAGAGGTAGCCCATCGCTTCAAGGATTTAGGGGAAGAAAACTTTAAAGCGCTGGTTCTGATTGCGTTTGCCCAATATTTGCAG
>2
TCAGAAGTAGCACATCGTTTCAAAGACCTTGGTGAAGAAAATTTTAAAGCTCTTGTTCTGATAGCTTTTGCTCAGTATCTGCAGCAGTGCCCGTTTGAAGATCATGTGAAACTGGTTAATGAAGTGACCGAATTTGCGAAAACCTGCGTCGCTGATGAGAGTGCCGAA
>3
CAGTGCCCGTTCGAAGACCATGTAAAATTAGTGAACGAAGTTACAGAGTTTGCCAAAACCTGCGTCGCGGATGAGTCAGCAGAAAATTGTGACAAGAGCCTGCATACCCTGTTTGGCGACAAACTGTGTACAGTCGCAACATTGCGCGAAACATATGGTGAGATGGCA
>4
AATTGCGATAAGAGCCTGCATACGTTGTTTGGTGACAAACTGTGCACCGTGGCAACGCTGCGTGAGACGTATGGGGAAATGGCTGATTGTTGTGCGAAACAAGAACCGGAGCGTAATGAGTGCTTTTTACAACATAAAGACGACAATCCGAATCTTCCGCGTTTAGTA
>5
GACTGTTGTGCTAAACAGGAACCGGAACGTAACGAATGCTTTCTGCAACATAAGGATGATAACCCGAACCTTCCGCGTTTAGTTCGTCCGGAAGTTGATGTGATGTGCACTGCATTCCATGACAATGAAGAGACTTTCCTGAAAAAGTATCTGTACGAAATAGCGCGC
```



######	Create Our Bowtie reference ( Elledge method )

```
bowtie-build our_vir3.fna our_vir3
```


Following their method, keep only first 50bp of reference.
I'm not sure I like this idea, but I'll give it a test.

```
cut -c1-50 our_vir3.fna > our_vir3.50.fna
bowtie-build our_vir3.50.fna our_vir3.50
```




######	Create Our Virus Score thresholds ( Elledge method )


Temporarily modify `scripts/elledge_calc_scores_nofilter.py` to reset the assigned peptides after each virus


```
echo "id,test_sample" > our_max_virus_scores.csv
for i in $( seq $( tail -n +2 our_vir3.csv | wc -l ) ) ; do
 echo "${i},True" >> our_max_virus_scores.csv
done

./scripts/elledge_calc_scores_nofilter.py max_virus_scores.csv our_vir3.csv Species 7 > our_max_virus_scores.virus_scores.csv
```


Our dataset is rather simple

```
echo "virus,count" > peptides_in_our_vir3.csv 
cat our_vir3.csv | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")";}(NR>1){print $3}' | sort | uniq -c | sed -e 's/^ *//' -e 's/ /,/' | awk -F, '{print $2","$1}' | sort -t, -k1,1 >> peptides_in_our_vir3.csv
```


And the thresholds are all simply 1.
```
echo "Virus,threshold" > our_vir3.virus_thresholds.csv
awk 'BEGIN{FS=OFS=","}(NR>1){m=$2/250;m=(m<1)?1:m;print $1,m}' our_max_virus_scores.virus_scores.csv > our_vir3.virus_thresholds.csv
```











###	Post Sequencing Analysis





####	Post (Larman method)

Scripts based on Larman paper's scripts or instructions.

Still just Pseudocode

#####	[Individual sample](scripts/larman_process_sample.bash)

PENDING

#####	[Aggregation of all samples](scripts/larman_aggregate_samples.bash)

PENDING







####	Post (Elledge method)

Scripts based on Elledge paper's scripts or instructions.







In “script.align.sh”, “bowtie -3 25” trims 25 nucleotides off the 3’ end of each sequencing read. This is done if sequencing reads are 75 nucleotides in length. The reference file only includes the first 50 nucleotides of each member of the library, so the sequencing reads must be trimmed down to 50 nucleotides to align correctly to the reference.





#####	Individual

######	[Align and extract counts for each FASTQ file](scripts/elledge_process_sample.bash)

* INPUT (FASTQ FILE)


```
for f in Elledge/fastq_files/*fastq.gz ; do
  s=$( echo $(basename $f .fastq.gz)| cut -d_ -f3 )
  elledge_process_sample.bash --sample_id ${s} --index ~/github/ucsffrancislab/PhIP-Seq/Elledge/vir3 ${f} 
done
```


#####	Aggregation

######	[sum all counts files for each sample](scripts/sum_counts_files.py)

* --output COMBINED_COUNTS_CSV_FILE
* List or glob of all counts files for the sample


```
for f in Elledge/fastq_files/LIB0*_L001_R1.count.csv.gz ; do
  i=${f/L001/L00*}
  o=${f/L001/L001_2_3_4}
  sum_counts_files.py --output ${o} ${i}
done
```



######	[Merge all combined counts files](scripts/merge_all_combined_counts_files.py)

* --output MERGED_COUNTS_CSV_FILE
* List or glob of all combined counts files



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



######	[Calculate Z-scores](scripts/elledge_Zscore_analysis.R)

* Merged combined counts csv file


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




######	[Booleanize Zscores at threshold of 3.5 for each technical replicate](scripts/booleanize_Zscore_replicates.py)

1. List of all technical replicates
2. Z-scores file



*QUESTION : Why is a threshold of 3.5 is used?*


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

######	[Merge booleanized Z-scores](scripts/merge_booleanized_replicates.py)

* --output MERGED_COMBINED_COUNTS_ZSCORE_BOOLEANIZED_REPLICATES_CSV_FILE
* List or glob of all replicates zscores combined counts files


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

This output file doesn't seem to be used elsewhere or usefull.


######	[Calculate virus scores](scripts/elledge_calc_scores_nofilter.py)

* hits - merged replicate zscores ( merged.combined.count.Zscores.S148.csv )
* oligo_metadata - Elledge/VIR3_clean.csv.gz ( we'll need to make one of these for our data )
* level - Species
* epitope_len - 7


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


######	Determining virus seropositivity

A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file “VirScan_viral_thresholds” contains the thresholds for each virus (Supplementary materials).
*NOTE : Public epitope annotations are available upon request.*

Based on the above line from the paper, I think that the following is accurate.



The virus score is the number of "novel" peptides for a given virus Species marked TRUE.
Novel determined by the longest match being less than the passed epitope_len (7)
Larger epitope_len would create higher virus scores.

This means that each virus score is dependent on the number of "novel" peptides.


This is why they use a separate threshold for each virus. 
Not sure where these treshold came from.

*NOTE : This threshold file only has thresholds for 206 viruses while the reference contains more than 440.*


*QUESTION : Why only have of the viruses in thresholds? Other irrelevant?*

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



######	Max virus scores

The virus score is the number of "novel" peptides for a given virus Species marked TRUE.
Novel determined by the longest match being less than the passed epitope_len (7)
Larger epitope_len would create higher virus scores.

This means that each virus score is dependent on the number of "novel" peptides.

``` 
tail -n 2 Elledge/vir3.fasta

>128287
ATTATGACAAGTTCAAAATTTGGCGGGGTCAATGTTTGGAATCGCTACTA
```

```
echo "id,test" > max_virus_scores.csv
for i in $( seq 128287 ) ; do
 echo "${i},True" >> max_virus_scores.csv
done

elledge_calc_scores_nofilter.py max_virus_scores.csv Elledge/VIR3_clean.csv.gz Species 7 > max_virus_scores.virus_scores.csv
```

58 viruses have a 0 max virus score? Something is off here.
```
awk -F, '($2==0)' max_virus_scores.virus_scores.csv  | wc -l
58

grep "^Adeno-associated virus," max_virus_scores.virus_scores.csv
Adeno-associated virus,0
```


```
echo "id,test" > test_virus_scores.csv
zgrep ",Adeno-associated virus," Elledge/VIR3_clean.csv.gz | awk -F, '{print $1",True"}' >> test_virus_scores.csv
elledge_calc_scores_nofilter.py test_virus_scores.csv Elledge/VIR3_clean.csv.gz Species 7 > test_virus_scores.virus_scores.csv

cat test_virus_scores.virus_scores.csv
Species,test
Adeno-associated virus,13
```
Hmm. Confused.


```
paste -d, original_max_virus_scores.virus_scores.csv <( cut -d, -f2 modified_max_virus_scores.virus_scores.csv ) > merged_max_virus_scores.virus_scores.csv

join --header -t, modified_max_virus_scores.virus_scores.csv Elledge/VirScan_viral_thresholds.csv > modified_max_virus_scores.virus_scores_with_thresholds.csv 
```

After discussion, this appears to be an intended bias as it is processed by the number of Species hits.











Additionally, I think that the order of comparison of the peptide can impact the virus score as well.
Thinking that one peptide "matches" several.
For example ...
```
id,peptide
1,ABCDEFGHIJKLMNOPQRSTUVWXYZ
2,ABCDEFGABCDEFGBCDEFGABCDEF
3,HIJKLMNHIJKLMNHIJKLMNHIJKL
4,OPQRSTUOPQRSTUOPQRSTUOPQRS
5,TUVWXYZTUVWXYZTUVWXYZTUVWX
```
If we found all 4, in this order we get a virus score of 1.
If we move 1 to last, we would get a virus score of 4.



```

echo "virus,count" > peptides_in_vir3_clean.csv 
zcat Elledge/VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")";}(NR>1){print $12}' | sort | uniq -c | sed -e 's/^ *//' -e 's/ /,/' | awk -F, '{print $2","$1}' | sort -t, -k1,1 >> peptides_in_vir3_clean.csv &

sort -t, -k1,1 modified_max_virus_scores.virus_scores_with_thresholds.csv > modified_max_virus_scores.virus_scores_with_thresholds.sorted.csv 

join --nocheck --header -t, peptides_in_vir3_clean.csv modified_max_virus_scores.virus_scores_with_thresholds.csv | awk 'BEGIN{FS=OFS=","}(NR==1){print $0,"mythreshold"}(NR>1){m=$3/250;m=(m<1)?1:m;print $0,m}' > modified_max_virus_scores.virus_scores_with_thresholds.peptides.csv

```


Generate a threshold file for all in the dataset.
```
echo "Virus,threshold" > modified_max_virus_scores.thresholds.csv
awk 'BEGIN{FS=OFS=","}(NR>1){m=$2/250;m=(m<1)?1:m;print $1,m}' modified_max_virus_scores.virus_scores.csv | sort -t, -k1,1 >> modified_max_virus_scores.thresholds.csv
```

The provided versus my thresholds shows that there are about a dozen additional. Perhaps they were ignoring them on purpose.
```
join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3)' | wc -l
join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv modified_max_virus_scores.thresholds.csv | awk -F, '($2>$3)' | wc -l

join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3)' | wc -l
join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv modified_max_virus_scores.thresholds.csv | awk -F, '($2>$3)' | wc -l
```

```
sdiff -Ws <( join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3){print $1,$2}' ) <( join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S148.virus_scores.csv modified_max_virus_scores.thresholds.csv | awk -F, '($2>$3){print $1,$2}' )
							      >	Bos taurus (Bovine) 2
							      >	Chikungunya virus (CHIKV) 2
							      >	Lassa mammarenavirus 4
							      >	Lipomyces starkeyi (Oleaginous yeast) 3
							      >	Middle East respiratory syndrome coronavirus 4
							      >	Pegivirus A 6
							      >	Saimiriine herpesvirus 2 (SaHV-2) (Herpesvirus saimiri) 3
							      >	Staphylococcus aureus 2
							      >	Streptococcus pneumoniae 11
							      >	Zika virus (strain Mr 766) (ZIKV) 2
							      >	deleted (mostly T. cruzi) 10
```

```
sdiff -Ws <( join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv Elledge/VirScan_viral_thresholds.csv | awk -F, '($2>$3){print $1,$2}' ) <( join --header -t, Elledge/fastq_files/merged.combined.count.Zscores.S150.virus_scores.csv modified_max_virus_scores.thresholds.csv | awk -F, '($2>$3){print $1,$2}' )
							      >	Bos taurus (Bovine) 3
							      >	Chikungunya virus (CHIKV) 3
							      >	Cosavirus A 2
							      >	Eastern equine encephalitis virus (EEEV) (Eastern equine ence
							      >	Lassa mammarenavirus 3
							      >	Lipomyces starkeyi (Oleaginous yeast) 2
							      >	Middle East respiratory syndrome coronavirus 3
							      >	Pegivirus A 5
							      >	Saimiriine herpesvirus 2 (SaHV-2) (Herpesvirus saimiri) 4
							      >	Staphylococcus aureus 6
							      >	Streptococcus pneumoniae 10
							      >	Zika virus (strain Mr 766) (ZIKV) 4
							      >	deleted (mostly T. cruzi) 7
```





##	Misc


Make a reference from the entire oligo instead of just the first 50bp.

Trim the adapters.

They aren't always the same.

```
zcat VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")"}(NR>1){o=substr($18,0,16);print o}' | sort | uniq -c
  70798 aGGAATTCCGCTGCGT
  37805 aTGAATTCGGAGCGGT

   5885 GGAATTCCGCTGCGTA
   4042 GGAATTCCGCTGCGTC
   6168 GGAATTCCGCTGCGTG
   3559 GGAATTCCGCTGCGTT

zcat VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")"}(NR>1){o=substr($18,184);print o}' | sort | uniq -c
  19654 CAGGGAAGAGCTCGA

   5961 aCACTGCACTCGAGACa
  11535 aCAGGgaagagctcgaa
  12958 cCACTGCACTCGAGACa
  23869 cCAGGgaagagctcgaa
   7386 gCACTGCACTCGAGACa
  13852 gCAGGgaagagctcgaa
  11500 tCACTGCACTCGAGACa
  21542 tCAGGgaagagctcgaa
```



```
zcat VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")"}(NR>1){o=substr($18,17,168);print ">"$17;print o}' | gzip > VIR3_clean.fna.gz
```

Not sure why, but some of these are repeated, hence the uniq/sort/uniq

```
zcat VIR3_clean.csv.gz | awk 'BEGIN{FPAT="([^,]*)|(\"[^\"]+\")"}(NR>1){o=substr($18,17,168);print $17","o}' | uniq | sort -t, -k1n,1 | uniq | awk -F, '{print ">"$1;print $2}' | gzip > VIR3_clean.uniq.fna.gz
```

```
bowtie-build VIR3_clean.uniq.fna.gz VIR3_clean
bowtie2-build VIR3_clean.uniq.fna.gz VIR3_clean
```








##	References, Papers & Resources

###	[PhIP-Seq: PHAGE DISPLAY IMMUNOPRECIPITATION SEQUENCING](https://cdi.bio/phip-seq/)

###	[Elledge 2022 - VirScan: High-throughput Profiling of Antiviral Antibody Epitopes](https://en.bio-protocol.org/en/bpdetail?id=4464&type=0)

###	[Elledge 2020 - Viral epitope profiling of COVID-19 patients reveals cross-reactivity and correlates of severity](https://www.science.org/doi/10.1126/science.abd4250)

https://doi.org/10.21769/BioProtoc.4464

https://os.bio-protocol.org/attached/file/20220629/Supplementary%20materials.docx

1) Videos: https://www.dropbox.com/sh/enlvqmsl1971bv3/AABfNDV21uO_bboLgWy3XBi4a?dl=0
2) Analysis scripts and input files: https://www.dropbox.com/sh/qvo1t75sgsq7fi8/AAAY-LQEQDrxV6wWF6OJDHPWa?dl=0
3) Example data: https://www.dropbox.com/sh/7hnqnx4yiskabgo/AAAeRE70EScRIzusqGX2HndUa?dl=0
4) Liquid handling robot protocol: https://www.dropbox.com/sh/guguiqgrmviqk9u/AABPoVam-Xsgg82su-IibDS_a?dl=0

This repo includes these analysis script, input files and example data.


###	[Larman 2018 - PhIP-Seq characterization of serum antibodies using oligonucleotide-encoded peptidomes](https://www.nature.com/articles/s41596-018-0025-6)

https://doi.org/10.1038/s41596-018-0025-6

https://github.com/lasersonlab/phip-stat

https://github.com/lasersonlab/pepsyn

The Larman paper did not include any specific code or data that isn't already easily accessible.


