#!/usr/bin/env python3

import sys
import pandas as pd


print(sys.argv)
print(sys.argv[1:3])
print(len(sys.argv))
print(sys.argv[1:len(sys.argv)-1])

threshold=3.5

columns=sys.argv[1:len(sys.argv)-1]
infile=sys.argv[len(sys.argv)-1]

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
outfile = outfile + "." + first_column + ".csv"

if( len(columns) > 0 ):
	out=df[[first_column]]>threshold
	print(out)

	if( len(columns) > 1 ):
		for column in columns[1:]:
			out=( ( out[first_column] ) & ( df[column]>threshold ) )
			print(outfile)
			pd.DataFrame( data=out, columns=[first_column]).to_csv(outfile)

