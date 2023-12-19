#!/usr/bin/env Rscript


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

datZ = virScanR::vs.convert(
		data = counts, 
		paramsList = convert_params1
	)


output_base = paste( inbase, "hits_combined_vir3_3.5_cutoff", "csv", sep=".")

mmR::mm.fastwrite(datZ$out, path = output_base )

