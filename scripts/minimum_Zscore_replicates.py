#!/usr/bin/env python3

import os
import sys
import pandas as pd

#print(sys.argv)
#print(sys.argv[1:3])
#print(len(sys.argv))
#print(sys.argv[1:len(sys.argv)-1])

#columns=sys.argv[1:len(sys.argv)-1]
#infile=sys.argv[len(sys.argv)-1]

# include standard modules
import argparse

# initiate the parser
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

#parser.add_argument('files', nargs='*', help='files help') #	list of samples name not files. really columns
parser.add_argument('columns', nargs='*', help='columns help. the column in the matrix to be treated as technical replicates')
parser.add_argument('-V','--version', action='version', version='%(prog)s 1.1')
parser.add_argument('-s', '--sample', nargs=1, type=str,
	help='output sample name to use %(prog)s (default: %(default)s)',required=True)
parser.add_argument('-m', '--matrix', nargs=1, type=str,
	help='zscore matrix filename to use %(prog)s (default: %(default)s)',required=True)
parser.add_argument('-o', '--output', nargs=1, type=str, default=['output.csv'],
	help='output csv filename %(prog)s (default: %(default)s)',required=True)

# read arguments from the command line
args = parser.parse_args()



if( len(args.columns) > 0 ):

	df = pd.read_csv(args.matrix[0],index_col="id")

	df[args.columns].min(axis=1).rename(args.sample[0]).to_csv(args.output[0])

else:
	print("No columns requested")


