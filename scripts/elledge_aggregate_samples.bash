#!/usr/bin/env bash

#	From section "Data analysis" B and C of "VirScan: High-throughput Profiling of Antiviral Antibody Epitopes"



#	Merge individual csv files into single table

#	merge into a table like ...

#	id,S148,S150,S154,S156,input
#	1,0,4,1,0,4
#	2,39,55,58,41,266
#	3,0,4,7,0,18
#	4,6,15,10,9,48
#	5,86,71,85,93,346

#	Will they always have the same number of sorted? Probably but don't count on it.






#	
#	B. Calculate Z-scores
#	
#	Note: To perform the Z-score analysis, count.combined files are merged into a table, and columns corresponding with no-serum controls are summed in a column called "input".
#	
#	1. Edit the R script "Zscore_analysis.R" to include the path to the count.combined table file and the desired
#	path to the output file, then run the script (Supplementary materials). The packages "mmR_0.1.0" and "virScanR_0.1.0.9000" are required (Supplementary materials).
#	Note: The file "Zscores_vir3" contains the results after this step (Supplementary materials).
#	2. A Z-score of at least 3.5 in both technical replicates of a sample is required to call a peptide a "hit".
#	Note: The file "hits_combined_vir3_3.5_cutoff" contains the results after this step (Supplementary materials).


#/francislab/data1/refs/refseq/phipSeq-20221116/Analysis_scripts_and_input_files/Zscore_analysis.R ${a_different_csv}






#	Convert this table into a boolean table above the threshold





#	C. Calculate virus scores
#	
#	1. Create a directory called "hits". In this directory should be .csv files for each sample with "True" or "False" values for each peptide ID, depending on whether the peptide scored as a hit (Z-score > 3.5) in both technical replicates of a sample or not. These files may be created by splitting each column of the "hits_combined_vir3_3.5_cutoff" file into a separate files (Supplementary materials).
#	2. Generate virus scores files using the following code:
#	Note: The "VIR3_clean" file provides the annotations for the oligos" (Supplementary materials). There are 115,753 oligos in the Vir3 library. Some protein fragments are identical in different viruses, and in these case there are multiple rows in the "VIR3_clean" file that correspond to a single oligo. To identify the viral source of a given peptide, look for the row(s) in the VIR3_clean file with the "id" value of the given peptide.











#	```
#	for i in hits/*.csv.gz; do python calc_scores_nofilter.py $i VIR3_clean.csv.gz Species 7 >virus_scores_$i; done
#	```
#	
#	D. Determining virus seropositivity
#	
#	1. A sample is determined to be seropositive for a virus if the virus_score > VirScan_viral_threshold and if at least one public epitope from that virus scores as a hit. The file "VirScan_viral_thresholds" contains the thresholds for each virus (Supplementary materials).
#	Note: Public epitope annotations are available upon request.
#	
#	



