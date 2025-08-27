#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy






# include standard modules
import argparse

# initiate the parser
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))

parser.add_argument('files', nargs='*', help='files help')
parser.add_argument('-V','--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-i', '--input', nargs=1, type=str, default=['Counts.csv'], help='input count csv filename to %(prog)s (default: %(default)s)')
#parser.add_argument('-o', '--output', nargs=1, type=str, default=['output.csv'], help='output zscore csv filename to %(prog)s (default: %(default)s)')

# read arguments from the command line
args = parser.parse_args()







print("Reading csv")
#df = pd.read_csv( args.input[0],index_col='id')
#df = pd.read_csv("ALL.csv",
df = pd.read_csv(args.input[0],
	header=[0,1,2,3,4,5,6,7,8],
	index_col=[0,1])
#	index_col=[0,1,2,3,4,5,6])

#df.columns = df.columns.set_names(["sample","subject","type","study","group","age","sex","plate"])
df.columns = df.columns.set_names(["sample","subject","bampath","type","study","group","age","sex","plate"])

#df.index = df.index.set_names(["id","species","protein","entry version","sequence version","start","end"])
df.index = df.index.set_names(["id","species"])

#df=df.droplevel([3,4,6])

df.shape







#	Why? Just Blank50_1?

#	Blank50_1 is on 20241204-Illumina-PhIP plate13 ( I should remove it from the manifest and not put this here )
#	Blank03_2 is on 20241204-Illumina-PhIP plate13 

#DO NOT USE Blank24dup plate6

#if( df.columns.get_level_values('sample').isin(['Blank50_1']).any() ):
#	df = df.drop('Blank50_1', level='sample', axis=1)








df.shape

# Create a function to replace the index values
def replace_unnameds(value):
	if value.startswith("Unnamed"):
		return("_")
	else:
		return value

df=df.rename(columns=replace_unnameds)

#df=df.sort_index(axis=1,level=[2,3,4,1,0])
df=df.sort_index(axis=1,level=[0,1,3,4,5])
df=df.fillna(0)

print(df.head())


samples=df.columns.get_level_values('sample').unique().to_numpy()
df_corr = pd.DataFrame(
	columns=samples,
	index=samples)

for sample1 in samples:
	for sample2 in samples:
		corr, p1 = scipy.stats.pearsonr(
			df.loc[:,df.columns.get_level_values('sample') == sample1],
			df.loc[:,df.columns.get_level_values('sample') == sample2]
		)
		df_corr.at[sample1,sample2] = corr

print(df_corr.head())
df_corr.to_csv( os.path.splitext( args.input[0] )[0]+".pearson.csv", sep=',', index=True)




print("Sums")
sums=df.sum()
print(sums.head())

print("Normalized CPM")
n=df/sums*1000000
print(n.head())

n=n.fillna(0)
#	use the "reset_index" technique otherwise, the multi-index ends up on a line by itself.
n.reset_index().to_csv( os.path.splitext( args.input[0] )[0]+".normalized.csv", sep=',', index=False)



###	print("Plate 1 Inputs")
###	print(n.loc[:,(n.columns.get_level_values('plate') == '1') & (n.columns.get_level_values('type') == 'input')])
###	
###	#print("Plate 1 Inputs Column Sums")
###	#print(n.loc[:,(n.columns.get_level_values('plate') == '1') & (n.columns.get_level_values('type') == 'input')].sum())
###	
###	print("Plate 1 Inputs Row Medians")
###	print(n.loc[:,(n.columns.get_level_values('plate') == '1') & (n.columns.get_level_values('type') == 'input')].median(axis=1))
###	
###	#print("Plate 1 Inputs Row Means")
###	#print(n.loc[:,(n.columns.get_level_values('plate') == '1') & (n.columns.get_level_values('type') == 'input')].mean(axis=1))
###	
###	s=n.sub(n.loc[:,(n.columns.get_level_values('plate') == '1') & (n.columns.get_level_values('type') == 'input')].median(axis=1),axis=0)
###	print(s.head())
###	
###	
###	print("Sum of abs diffs")
###	print(( df.loc[:,df.columns.get_level_values('sample') == 'CSE01_1'].droplevel('sample',axis=1) - df.loc[:,df.columns.get_level_values('sample') == 'CSE01_2'].droplevel('sample',axis=1)).abs().sum())
###	print(( n.loc[:,n.columns.get_level_values('sample') == 'CSE01_1'].droplevel('sample',axis=1) - n.loc[:,n.columns.get_level_values('sample') == 'CSE01_2'].droplevel('sample',axis=1)).abs().sum())
###	print(( s.loc[:,s.columns.get_level_values('sample') == 'CSE01_1'].droplevel('sample',axis=1) - s.loc[:,s.columns.get_level_values('sample') == 'CSE01_2'].droplevel('sample',axis=1)).abs().sum())
###	






print(n.columns.get_level_values('plate').unique())

ss=[]
for plate in n.columns.get_level_values('plate').unique():
	a = n.loc[:,(n.columns.get_level_values('plate') == plate)]
	ss.append(
		a.sub(
			a.loc[:,(a.columns.get_level_values('type') == 'input')].median(axis=1),
				axis=0)
	)

s=pd.concat(ss,axis=1)
s=s.sort_index(axis=1)
s=s.fillna(0)

print( s.head() )

#s.to_csv('s5.csv', sep=',', index=True, header=True)
#s.to_csv( os.path.splitext( args.input[0] )[0]+".normalized.subtracted.csv", sep=',', index=True, header=True)
#	use the "reset_index" technique otherwise, the multi-index ends up on a line by itself.
s.reset_index().to_csv( os.path.splitext( args.input[0] )[0]+".normalized.subtracted.csv", sep=',', index=False)

samples=s.columns.get_level_values('sample').unique().to_numpy()
s_corr = pd.DataFrame( columns=samples, index=samples )
for sample1 in samples:
	for sample2 in samples:
		corr, p1 = scipy.stats.pearsonr(
			s.loc[:,s.columns.get_level_values('sample') == sample1],
			s.loc[:,s.columns.get_level_values('sample') == sample2]
		)
		s_corr.at[sample1,sample2] = corr
	
s_corr=s_corr.fillna(0)
#print(s_corr)
#s_corr.to_csv('s_corr5.csv', sep=',', index=True, header=True)
s_corr.to_csv( os.path.splitext( args.input[0] )[0]+".normalized.subtracted.pearson.csv", sep=',', index=True)







#	zf = pd.read_csv("zscores.csv",
#		header=[0,1,2,3,4,5,6,7],
#		index_col=[0])
#	zf.columns = zf.columns.set_names(["sample","subject","type","study","group","age","sex","plate"])
#	zf.index = zf.index.set_names(["id"])
#	#zf=zf.sort_index(axis=1,level=[2,3,4,1,0])
#	zf=zf.sort_index(axis=1,level=[0,1,2,3,4])
#	zf=zf.fillna(0)
#	
#	samples=zf.columns.get_level_values('sample').unique().to_numpy()
#	z_corr = pd.DataFrame( columns=samples, index=samples)
#	
#	print("Correlating each sample")
#	for sample1 in samples:
#		for sample2 in samples:
#			corr, p1 = scipy.stats.pearsonr(
#				zf.loc[:,zf.columns.get_level_values('sample') == sample1],
#				zf.loc[:,zf.columns.get_level_values('sample') == sample2]
#			)
#			z_corr.at[sample1,sample2] = corr
#	
#	print("Writing")
#	z_corr=z_corr.fillna(0)
#	print(z_corr)
#	z_corr.to_csv('z_corr4.csv', sep=',', index=True, header=True)
#	
#	
#	
#	
#	
#	
#	#	#	#	#z = pd.read_csv("../out.plate1/All.count.Zscores.csv",
#	#	#	#	z = pd.read_csv("../out.all/All.count.Zscores.csv",
#	#	#	#		header=[0],
#	#	#	#		index_col=[0])
#	#	#	#	z.columns = z.columns.set_names(["sample"])
#	#	#	#	z.index = z.index.set_names(["sample"])
#	#	#	#	z.fillna(0,inplace=True)
#	#	#	#	
#	#	#	#	z_corr = pd.DataFrame(
#	#	#	#		columns=z.columns.get_level_values('sample').to_numpy().unique(),
#	#	#	#		index=z.columns.get_level_values('sample').to_numpy().unique())
#	#	#	#	
#	#	#	#	for sample1 in z.columns.get_level_values('sample').to_numpy().unique():
#	#	#	#		for sample2 in z.columns.get_level_values('sample').to_numpy().unique():
#	#	#	#			#corr = np.correlate(
#	#	#	#			corr, p1 = scipy.stats.pearsonr(
#	#	#	#				z.loc[:,z.columns.get_level_values('sample') == sample1],
#	#	#	#				z.loc[:,z.columns.get_level_values('sample') == sample2]
#	#	#	#			)
#	#	#	#			z_corr.at[sample1,sample2] = corr
#	#	#	#	
#	#	#	#	print(z_corr)
#	#	#	#	z_corr.to_csv('z_corr.csv', sep=',', index=True, header=True)
#	#	
#	#	
#	#	
#	#	#		print(n.columns.get_level_values('plate').unique())
#	#	#		
#	#	#		print("Reading Zscores files")
#	#	#		zz=[]
#	#	#		for plate in n.columns.get_level_values('plate').unique():
#	#	#			#print("Reading ../out.plate"+plate+"/Zscores.t.csv")
#	#	#			#a = pd.read_csv("../out.plate"+plate+"/Zscores.t.csv",
#	#	#			#	header=[0,1,2],
#	#	#			#	index_col=[0,1])
#	#	#			print("Reading ../out.plate"+plate+"/All.count.Zscores.csv")
#	#	#			a = pd.read_csv("../out.plate"+plate+"/All.count.Zscores.csv",
#	#	#				header=[0],
#	#	#				index_col=[0])
#	#	#			print(a.head())
#	#	#			#a.loc[('plate',''),:]=plate
#	#	#			a.loc['plate',:]=plate
#	#	#			a.columns = pd.MultiIndex.from_arrays([
#	#	#				a.columns.get_level_values(0).to_numpy(),
#	#	#				a.loc['plate'].to_numpy() ])
#	#	#				#a.loc[('plate','')].to_numpy() ])
#	#	#				#a.columns.get_level_values(1).to_numpy(),
#	#	#				#a.columns.get_level_values(2).to_numpy(),
#	#	#			#a.drop(('plate',''),inplace=True)
#	#	#			a.drop('plate',inplace=True)
#	#	#			#a.columns = a.columns.set_names(["subject","type","sample","plate"])
#	#	#			#a.index = a.index.set_names(["id","species"])
#	#	#			a.columns = a.columns.set_names(["sample","plate"])
#	#	#			a.index = a.index.set_names(["id"])
#	#	#			zz.append(a)
#	#	#		
#	#	#		print("Concatting")
#	#	#		
#	#	#		z=pd.concat(zz,axis=1)
#	#	#		
#	#	#		print(z.columns)
#	#	#		print(z.columns.get_level_values('sample').unique().to_numpy())
#	#	#		
#	#	#		#z=z.sort_index(axis=1,level=[3,1,0,2])	#	sort columns plate, type, subject, sample
#	#	#		#z=z.sort_index(axis=1,level=[1,0,2])	#	sort columns type, subject, sample
#	#	#		z=z.sort_index(axis=1,level=[0])	#	sort columns sample
#	#	#		z=z.fillna(0)
#	#	#		
#	#	#		print("Writing")
#	#	#		print(z.head())
#	#	#		z.to_csv('z.csv', sep=',', index=True, header=True)
#	#	#		
#	#	#		samples=z.columns.get_level_values('sample').unique().to_numpy()
#	#	#		z_corr = pd.DataFrame(
#	#	#			columns=samples,
#	#	#			index=samples)
#	#	#		
#	#	#		print(z.columns.get_level_values('sample'))
#	#	#		
#	#	#		print("Correlating each sample")
#	#	#		for sample1 in samples:
#	#	#			for sample2 in samples:
#	#	#				#print( sample1, " - ", sample2 )
#	#	#				corr, p1 = scipy.stats.pearsonr(
#	#	#					z.loc[:,z.columns.get_level_values('sample') == sample1],
#	#	#					z.loc[:,z.columns.get_level_values('sample') == sample2]
#	#	#				)
#	#	#				z_corr.at[sample1,sample2] = corr
#	#	#			
#	#	#		print("Writing")
#	#	#		print(z_corr)
#	#	#		z_corr.to_csv('z_corr.csv', sep=',', index=True, header=True)


