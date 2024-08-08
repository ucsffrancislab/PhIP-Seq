#!/usr/bin/env python3

import os
import sys
import pandas as pd

print(sys.argv)
print(sys.argv[1:3])
print(len(sys.argv))
print(sys.argv[1:len(sys.argv)-1])

threshold=3.5

#columns=sys.argv[1:len(sys.argv)-1]
#infile=sys.argv[len(sys.argv)-1]

# include standard modules
import argparse

# initiate the parser
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

#parser.add_argument('files', nargs='*', help='files help') #	list of samples name not files. really columns
parser.add_argument('columns', nargs='*', help='columns help. the column in the matrix to be treated as technical replicates')
parser.add_argument('-V','--version', action='version', version='%(prog)s 1.1')
#parser.add_argument('-o', '--output', nargs=1, type=str, default=['merged.csv.gz'], help='output csv filename %(prog)s (default: %(default)s)')
parser.add_argument('-s', '--sample', nargs=1, type=str, help='output sample name to use %(prog)s (default: %(default)s)')
parser.add_argument('-m', '--matrix', nargs=1, type=str, help='zscore matrix filename to use %(prog)s (default: %(default)s)')

# read arguments from the command line
args = parser.parse_args()

#sample=args.sample[0]
#print( "Using sample name: ", sample )

columns=args.columns
infile=args.matrix[0]


#from pathlib import Path
#p = Path(infile)
#extensions = "".join(p.suffixes)

df = pd.read_csv(infile,index_col="id")

#	not sure if necessary
#df.drop(["group","input"],axis="columns",inplace=True)


first_column=columns[0]

import re
#outfile = infile.replace(".gz$", "")
#outfile = infile.replace(".csv$", "")
outfile = re.sub('.gz$', '', infile)
outfile = re.sub('.csv$', '', infile)
outfile = outfile + "." + args.sample[0] + ".csv"

if( len(columns) > 0 ):
	out=df[[first_column]]>threshold
	print(out)

	if( len(columns) > 1 ):
		for column in columns[1:]:
			out=( ( out[first_column] ) & ( df[column]>threshold ) )
			print(outfile)
			pd.DataFrame( data=out, columns=[first_column]).to_csv(outfile)

