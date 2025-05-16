#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy

#	phip_seq_phagelib_normalize_counts.py -i out.plate1/Counts.csv




# include standard modules
import argparse

# initiate the parser
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

#parser.add_argument('files', nargs='*', help='files help')
parser.add_argument('-V','--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-i', '--input', nargs=1, type=str, default=['Counts.csv'], help='input count csv filename to %(prog)s (default: %(default)s)')
#parser.add_argument('-o', '--output', nargs=1, type=str, default=['output.csv'], help='output zscore csv filename to %(prog)s (default: %(default)s)')

# read arguments from the command line
args = parser.parse_args()







print("Reading csv")
df = pd.read_csv(args.input[0],
	header=[0,1,2,3,4,5,6,7,8],
	index_col=[0,1])

df.columns = df.columns.set_names(["sample","subject","bampath","type","study","group","age","sex","plate"])

df.index = df.index.set_names(["id","species"])

df.shape





#	Why? Just Blank50_1?

#	Blank50_1 is on 20241204-Illumina-PhIP plate13 ( I should remove it from the manifest and not put this here )

#DO NOT USE Blank24dup plate6

if( df.columns.get_level_values('sample').isin(['Blank50_1']).any() ):
	df = df.drop('Blank50_1', level='sample', axis=1)





df.shape

# Create a function to replace the index values
def replace_unnameds(value):
	if value.startswith("Unnamed"):
		return("_")
	else:
		return value

df=df.rename(columns=replace_unnameds)

df=df.sort_index(axis=1,level=[0,1,3,4,5])
df=df.fillna(0)

print(df.head())


samples=df.columns.get_level_values('sample').unique().to_numpy()
#	df_corr = pd.DataFrame(
#		columns=samples,
#		index=samples)
#	
#	for sample1 in samples:
#		for sample2 in samples:
#			corr, p1 = scipy.stats.pearsonr(
#				df.loc[:,df.columns.get_level_values('sample') == sample1],
#				df.loc[:,df.columns.get_level_values('sample') == sample2]
#			)
#			df_corr.at[sample1,sample2] = corr
#	
#	print(df_corr.head())
#	df_corr.to_csv( os.path.splitext( args.input[0] )[0]+".pearson.csv", sep=',', index=True)




print("Sums")
sums=df.sum()
print(sums.head())

print("Normalized CPM")
n=df/sums*1000000
print(n.head())

n=n.fillna(0)




#	use the "reset_index" technique otherwise, the multi-index ends up on a line by itself.
#n.reset_index().to_csv( os.path.splitext( args.input[0] )[0]+".normalized.csv", sep=',', index=False)



#	https://www.science.org/doi/10.1126/science.adc9498
#	We first mapped the sequencing reads to the reference library sequences using Bowtie (55) and counted the number of reads corresponding to each peptide in the input library and each sample “output”. For each sample, we normalized the read counts for each peptide by the total read counts for the sample. Then, we divided the normalized read counts of each peptide in the sample by the normalized read counts of each peptide in the input library to obtain an enrichment value. We averaged enrichment values for technical replicates of a sample.
#df[df.loc[:,(df.columns.get_level_values('type') == 'Phage Library')].median(axis='columns')>0]

print(n.columns.get_level_values('plate').unique())

ss=[]
for plate in n.columns.get_level_values('plate').unique():
	a = n.loc[:,(n.columns.get_level_values('plate') == plate)]
	a = a[a.loc[:,(a.columns.get_level_values('type') == 'Phage Library')].median(axis='columns')>0]
	ss.append(
		a.div(
			a.loc[:,(a.columns.get_level_values('type') == 'Phage Library')].median(axis=1),
				axis=0)
	)

s=pd.concat(ss,axis=1)
s=s.sort_index(axis=1)
s=s.fillna(0)

print( s.head() )

#	use the "reset_index" technique otherwise, the multi-index ends up on a line by itself.
s.reset_index().to_csv( os.path.splitext( args.input[0] )[0]+".normalized.normalized.csv", sep=',', index=False)



s.droplevel([2,4,5,6,7,8],axis='columns').reorder_levels(["subject","type","sample"],axis='columns').T.reset_index().to_csv( os.path.splitext( args.input[0] )[0]+".normalized.normalized.trim.csv", sep=',', index=False)








#	
#	samples=s.columns.get_level_values('sample').unique().to_numpy()
#	s_corr = pd.DataFrame( columns=samples, index=samples )
#	for sample1 in samples:
#		for sample2 in samples:
#			corr, p1 = scipy.stats.pearsonr(
#				s.loc[:,s.columns.get_level_values('sample') == sample1],
#				s.loc[:,s.columns.get_level_values('sample') == sample2]
#			)
#			s_corr.at[sample1,sample2] = corr
#		
#	s_corr=s_corr.fillna(0)
#	s_corr.to_csv( os.path.splitext( args.input[0] )[0]+".normalized.normalized.pearson.csv", sep=',', index=True)
#	

