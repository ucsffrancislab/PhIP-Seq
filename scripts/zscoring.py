#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats



# include standard modules
import argparse

# initiate the parser
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

parser.add_argument('files', nargs='*', help='files help')
parser.add_argument('-V','--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-i', '--input', nargs=1, type=str, default=['All.count.csv'], help='input count csv filename to %(prog)s (default: %(default)s)')
parser.add_argument('-o', '--output', nargs=1, type=str, default=['output.csv'], help='output zscore csv filename to %(prog)s (default: %(default)s)')

#parser.add_argument('--count_threshold', type=int, default=300, help='Tile count minimum threshold to %(prog)s (default: %(default)s)')

#parser.add_argument('-n', '--normal', nargs=1, type=str, required=True, help='normal matrix  filename to %(prog)s (default: %(default)s)')
#parser.add_argument('-t', '--tumor',  nargs=1, type=str, required=True, help='tumor matrix  filename to %(prog)s (default: %(default)s)')
#parser.add_argument('-s', '--select', nargs=1, type=str, required=True, help='[GENE,Symbol,RE] select list  filename to %(prog)s (default: %(default)s)')
#parser.add_argument('-o', '--output', nargs=1, type=str, default=['output.tsv'], help='output tsv filename to %(prog)s (default: %(default)s)')

# read arguments from the command line
args = parser.parse_args()







print("Reading csv")
df = pd.read_csv( args.input[0],index_col='id')

print("\nFirst 10x10")
print(df.iloc[:10,:10] )

print("\nLast 10x10")
print(df.iloc[-10:,-10:] )

samples = df.columns.to_list()
samples.remove("input")
print("\nSamples:")
print(samples)

print("\nExtracting input to list")
inputs=df['input'].to_list()
print("\nFirst 10 input")
print(inputs[:10])
#[38, 0, 1170, 197, 0]

print("\nConverting to array")
sorted_inputs=np.sort(inputs).tolist()
print(sorted_inputs[0:5])
#[0, 0, 0, 0, 0]
print(sorted_inputs[(len(sorted_inputs)-5):(len(sorted_inputs))])
#[100553, 101219, 124937, 134726, 274421]

print("\nCounting occurences and saving in dict")
counts=dict((keys, sorted_inputs.count(keys)) for keys in sorted_inputs)
#print(d)
#{0: 16728, 1: 6996, 2: 3479, 3: 2469, 4: 2036, 5: 1898, 6: 1684, 7: 1458, 
#...
#, 100553: 1, 101219: 1, 124937: 1, 134726: 1, 274421: 1}

keys=counts.keys()
print(type(keys))
#<class 'dict_keys'>
print(type(list(keys)[0]))
#<class 'int'>
print(list(keys)[0])
#0

print("\nSorting keys")
sorted_keys=sorted(list(keys),reverse=True)
print( sorted_keys[0:9] )
#	[274421, 134726, 124937, 101219, 100553, 99564, 65635, 64566, 63220]


#	cut -d, -f46 All.count.csv | tail -n +2 | sort -n | uniq -c | sed 's/^\s*//g' | sort -k2nr,2 | awk '{if(s>=300){print s,list;s=0;list="";}s+=$1;list=list","$2}END{print s,list}'



#	cut -d, -f46 All.count.csv | tail -n +2 | sort -n | uniq -c | sed 's/^\s*//g' | sort -k2nr,2 | head
#	1 274421
#	1 134726
#	1 124937
#	1 101219
#	1 100553
#	1 99564
#	1 65635
#	1 64566
#	1 63220
#	1 62560
#	
#	cut -d, -f46 All.count.csv | tail -n +2 | sort -n | uniq -c | sed 's/^\s*//g' | sort -k2nr,2 | tail
#	1353 9
#	1452 8
#	1458 7
#	1684 6
#	1898 5
#	2036 4
#	2469 3
#	3479 2
#	6996 1
#	16728 0

#	If we start at the big end, we run the risk of only having a few in the last bin
#	If we start at the low end, we can get a VERY broad spread of noise. Preferred.
#	Does that matter?
#	Should we consider the spread as well as the count?


print("\nBinning")
tile_count=0
bin=1
df['bin']=0
max_input=0
for key in sorted_keys:
	if max_input == 0:
		max_input=key

	input_diff = max_input - key
	if input_diff > 1500:
		count_threshold=50
	elif input_diff > 1000:
		count_threshold=75
	elif input_diff > 500:
		count_threshold=150
	else:
		count_threshold=300
			
	#	if( tile_count >= args.count_threshold ):
	if( tile_count >= count_threshold ):
		print("Tile count surpassed progressive count threshold:",count_threshold)
		bin+=1
		#print( tile_count )
		tile_count=0
		max_input=key

	tile_count+=counts[key]
	df.loc[df['input']==key,'bin']=bin


print("\nDone\n")

print("\nFirst 10x10")
print(df.iloc[:10,:10] )

print("\nLast 10x10")
print(df.iloc[-10:,-10:] )




#	"The formula for calculating a z-score is z = (x-μ)/σ, where x is the raw score, μ is the population mean, and σ is the population standard deviation."
#	Using the blanks, they rank and then cluster tiles with similar noise, use the middle ~80% or so of the tiles within a cluster to determine the cluster mean and standard deviation, and then for every tile in the cluster compute the z score from that specific set of parameters.
#	The z scores are computed across tiles for an individual sample. The other samples should have no influence on it.
#	The basics of the VirScan pipeline are:
#	Rank tiles based on the input column.
#	Bin neighboring tiles into groups of at least 300, where tiles with identical rank go into the same bin. These bins represent tiles with nearly identical background noise.
#	Creation of bins with uniform baseline epitope abundance: Next, bins were created by first binning epitopes with identical ranks and then, if the bin contained fewer than 300 epitopes, additional epitopes with adjacent ranks were added until each bin contained at minimum 300 epitopes. VirScan 2.0 contains approximately 115,000 epitopes and it was therefore common for > 300 epitopes to be tied for the same exact rank (fig. M1).



#	(The rest is done Per sample...)
#	For each bin, collect all of the tile values for sample i and collect the middle 90% of values.
#	From the middle 90%, calculate the mean and standard deviation.
#	For all tiles in the bin (including the left out upper and lower 5%'s, calculate Z score  = (x-mean)/s.d.


out = pd.DataFrame().reindex_like(df)
out['bin'] = df['bin']
out['input'] = df['input']


#def zscore(x,mean,stddev):
#	if x > 1:
#		return ( x - mean ) / stddev
#	else:
#		return ( 0 - mean ) / stddev
#out.loc[out['bin']==i,[sample]] = df.loc[df['bin']==i,[sample]].apply(



#/c4/home/gwendt/.local/bin/zscoring.py:191: RuntimeWarning: divide by zero encountered in scalar divide
#  x[x>1].apply( lambda y: ( y - mean ) / stddev )
#/c4/home/gwendt/.local/bin/zscoring.py:192: RuntimeWarning: invalid value encountered in scalar divide
#  x[x<=1].apply( lambda y: ( 0 - mean ) / stddev )

np.seterr(divide='ignore', invalid='ignore')


#def zscore(series,mean,stddev):	#	untested
#	return series.apply( lambda cell: ( ( cell if cell > 1 else 0 ) - mean ) / stddev )

#def zscore(series,mean,stddev):	#	works
#	z=x[x>1].apply( lambda y: ( y - mean ) / stddev )
#	z=x[x<=1].apply( lambda y: ( 0 - mean ) / stddev )
#	return z


#import pandas as pd
#
#df = pd.DataFrame({
#    'A': [1, 2, 3],
#    'B': [4, 5, 6]
#})
#
#def inner_apply(row):
#    return row['A'] * 2
#
#def outer_apply(row):
#    row['C'] = inner_apply(row)
#    return row
#
#df = df.apply(outer_apply, axis=1)
#
#print(df)




#lambda x: ( ( x if x > 1 else 0 ) - m ) / std )

#	return ( ( x if x > 1 else 0 ) - mean ) / stddev #	doesn't work!

#	return ( x - mean ) / stddev # will return +inf, -inf and NaN

#	return float('nan') if stddev == 0.0 else ( x - mean ) / stddev #	remove the +inf / -inf


#np.set_printoptions(threshold=np.inf)  # Set threshold to infinity
#pd.set_option('display.max_rows', None)
#pd.set_option('display.max_columns', None)


#Traceback (most recent call last):
#  File "/c4/home/gwendt/.local/bin/zscoring.py", line 218, in <module>
#    out.loc[out['bin']==i,[sample]] = df.loc[df['bin']==i,[sample]].apply(
#                                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  File "/c4/home/gwendt/.local/lib/python3.11/site-packages/pandas/core/frame.py", line 10374, in apply
#    return op.apply().__finalize__(self, method="apply")
#           ^^^^^^^^^^
#  File "/c4/home/gwendt/.local/lib/python3.11/site-packages/pandas/core/apply.py", line 916, in apply
#    return self.apply_standard()
#           ^^^^^^^^^^^^^^^^^^^^^
#  File "/c4/home/gwendt/.local/lib/python3.11/site-packages/pandas/core/apply.py", line 1063, in apply_standard
#    results, res_index = self.apply_series_generator()
#                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  File "/c4/home/gwendt/.local/lib/python3.11/site-packages/pandas/core/apply.py", line 1081, in apply_series_generator
#    results[i] = self.func(v, *self.args, **self.kwargs)
#                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  File "/c4/home/gwendt/.local/bin/zscoring.py", line 219, in <lambda>
#    lambda x: ( ( x if x > 1 else 0 ) - m ) / std )
#                       ^^^^^
#  File "/c4/home/gwendt/.local/lib/python3.11/site-packages/pandas/core/generic.py", line 1577, in __nonzero__
#    raise ValueError(
#ValueError: The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all().


for sample in samples:
#for sample in ['14627-01dup']:
	print("Sample:",sample)
	for i in range(1,bin+1):
	#for i in [146]:
		print("i:",i)

		#tmp=df.loc[df['bin']==i,[sample]]
		#mask=tmp>0
		#print(tmp[mask.loc[:,sample]])

		X = df.loc[df['bin']==i,[sample]].values.flatten()
		#	or .explode()
		#print(X)
		m = stats.trim_mean(X, 0.05)
		print("Trimmed mean:",m)
		std = stats.mstats.trimmed_std(X,0.05)
		print("Trimmed stddev:",std)

		#out.loc[out['bin']==i,[sample]] = df.loc[df['bin']==i,[sample]].apply(zscore,args=(m,std))

		#out.loc[out['bin']==i,[sample]] = df.loc[df['bin']==i,[sample]].apply(
		#	lambda x: ( ( x if x > 1 else 0 ) - m ) / std )

		out.loc[out['bin']==i,[sample]] = df.loc[df['bin']==i,[sample]].apply(
			lambda series: series.apply( 
				lambda cell: ( ( cell if ( cell > 1 ) else 0 ) - m ) / std ) )

		#df['C'] = df['A'].apply(lambda x: x * 2 if x > 1 else x)

		#print(out.loc[out['bin']==i,[sample]])


out.to_csv( args.output[0], index_label=['id'] )


