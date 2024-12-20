#!/usr/bin/env Rscript


#	take an input file like ...
#	head ~/testing/counts_vir3_Bio_Protocol.csv
#	id,S148,S150,S154,S156,input
#	1,0,4,1,0,4
#	2,39,55,58,41,266
#	3,0,4,7,0,18
#	4,6,15,10,9,48
#	5,86,71,85,93,346
#	6,0,4,5,4,12
#	7,54,59,44,77,5
#	8,5,11,6,9,31
#	9,1,0,3,0,7

#	and create a file like ... (where does "group" come from?)

#	head ~/testing/Zscores_vir3.csv
#	id,group,S148,S150,S154,S156,input
#	1,5,-0.71,2.99,0.21,-0.78,4
#	2,217,-0.36,0.24,1.24,-1.07,266
#	3,19,-1.32,0.35,1.99,-1.37,18
#	4,49,-0.31,1.74,0.59,-0.27,48
#	5,235,1.67,0.22,1.82,0.86,346
#	6,13,-1.02,1.04,1.77,0.75,12
#	7,6,52,44.15,40.07,52.93,5
#	8,32,0.16,1.88,0.47,0.8,31
#	9,8,-0.03,-0.9,1.59,-0.96,7



#	R CMD INSTALL mmR_0.1.0.tar.gz
#	R -e "install.packages('VGAM')"
#	R CMD INSTALL virScanR_0.1.0.9000.tar.gz


library(virScanR)
library(mmR)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
args[1] = normalizePath(args[1])

print("input")
print(args[1])


inbase=gsub(".gz","",args[1])
print("inbase")
print(inbase)
inbase=gsub(".csv","",inbase)
print(inbase)

#quit()

#Taking input= as a system command ('gunzip -c /francislab/data2/refs/refseq/phipSeq-20221116/testing/ENCFF878GFR.count.csv.gz') and a variable has been used in the expression passed to `input=`. Please use fread(cmd=...). There is a security concern if you are creating an app, and the app could have a malicious user, and the app is not running in a secure environment; e.g. the app is running as root. Please read item 5 in the NEWS file for v1.11.6 for more information and for the option to suppress this message.

counts <- mmR::mm.fastread(args[1])

#	Doesn't like gzipped csv ... but seems to work otherwise?
#Warning messages:
#1: In grepl("SQLite", a) :
#  unable to translate '<8/<a1>|e' to a wide string
#2: In grepl("SQLite", a) : input string 1 is invalid
#3: In grepl("SQLite", a) :
#  unable to translate '<8/<a1>|e' to a wide string
#4: In grepl("SQLite", a) : input string 1 is invalid
#5: In grepl("SQLite", a) :
#  unable to translate '<8/<a1>|e' to a wide string
#6: In grepl("SQLite", a) : input string 1 is invalid

head(counts)

propTrim = .04

#two other pertinent: NLP_nBinom", "NLP_pois"

convert_params1 <- virScanR::vs.set_params_convert(stat  = "Z_score",
		makeGroups = TRUE,
		makeGroupSize = 300,
		idCol = "id",
		groupTogetherIfGreaterThanGroupSize = TRUE,
		splitHighVarMeanGrps = TRUE,
		cols_to_remove = NULL,
		cols_to_not_evaluate = c("id","input"),
		propExtremesToRemove = propTrim,
		removeTail = "both",
		coresAcrossGroups = parallel::detectCores()-2,
		returnAs = "data.table",
		returnParams = TRUE)
print(convert_params1)

datZ = virScanR::vs.convert(
		data = counts,
		paramsList = convert_params1
	)

print(datZ)


output_base = paste( inbase, "Zscores", "csv", sep=".")

mmR::mm.fastwrite(datZ$out, path = output_base )

