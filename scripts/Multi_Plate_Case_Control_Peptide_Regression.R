#!/usr/bin/env Rscript

#	Uses a logistic regression framework to assess the association between the presence of each peptide and user specified case control status.
#	Uses a user specifed Z score threshold (i.e. Z >3.5) to call positivity
#	Takes as input a list of directories (e.g. plates), to use in the analysis, the directories should contain the Zscores.csv and the manifest* files.
#	Only peptides found in both files are included in the analysis
#	The logistic regression adjusts for age (continuous), sex (M/F factor), and plate (1..n plates, factor)
#	The output (to a user specified directory owd , is a list of peptides, ordered by ascending p-value. It includes peptide ID, species, frequency of peptide in the first group, frequency of the peptide in the second group, and the beta, standard error, and p-value from the multivariable logistic regression (NA if there is no variation in the presence of the peptide overall, e.g. all 0 or all 1).
#	The file will also output a logfile with similar naming convention, which includes the details of the plates used in the analysis, sample sizes, etc.

#	Expecting this precise format...
#	hc out.plate14/Counts.normalized.subtracted.trim.select-121314.csv 
#	y,x,id,1,10,100,1000,10000,10001,10002,10003,10004,10005,10006,10007,10008,10009,1001,10010,10011,10
#	subject,type,species,Papiine herpesvirus 2,Vaccinia virus,Human herpesvirus 3,Hepatitis B virus,Huma
#	024JCM,pemphigus serum,024JCM,0.0,0.0,9.16111056217155,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0
#	024JCM,pemphigus serum,024JCMdup,0.0,0.0,0.6039953081644461,0.6039953081644461,0.0,0.0,2.41598123265
#	074KBP,pemphigus serum,074KBP,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.118460014774857,0.0,0.0,0.0,0.0,0.0,0.0,
#	074KBP,pemphigus serum,074KBPdup,0.0,0.0,3.4158313535744114,0.0,6.831662707148823,1.1386104511914705
#	
#	hc out.plate14/Zscores.select-121314.csv 
#	y,x,id,1,10,100,1000,10000,10001,10002,10003,10004,10005,10006,10007,10008,10009,1001,10010,10011,10
#	subject,type,species,Papiine herpesvirus 2,Vaccinia virus,Human herpesvirus 3,Hepatitis B virus,Huma
#	024JCM,pemphigus serum,024JCM,-0.22559419625872576,-0.22559419625872576,2.1116455734686204,-0.225594
#	024JCM,pemphigus serum,024JCMdup,-0.3524648371111257,-0.3524648371111257,-0.5719156651485517,-0.3524
#	074KBP,pemphigus serum,074KBP,-0.3906228577599959,-0.3906228577599959,-0.6522456635017637,-0.3906228
#	074KBP,pemphigus serum,074KBPdup,-0.44436930651974654,-0.44436930651974654,-0.17707278654036396,-0.4


library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
scriptdir=dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(paste(scriptdir,'GenoLib.R',sep='/'))
parser <- ArgumentParser(description=scriptname)
parser$add_argument("-t", "--type", type="character", default="",
	help="limit type", metavar="type")
parser$add_argument("-a", "--group1", type="character", required=TRUE,
	help="first group to compare", metavar="group")
parser$add_argument("-b", "--group2", type="character", required=TRUE,
	help="second group to compare", metavar="group")
parser$add_argument("-z", "--zscore", type="double", default=3.5,
	help="zscore threshold", metavar="double")
parser$add_argument("-s", "--sex", type="character", default="",
	help="limit sex", metavar="sex")
parser$add_argument("-o", "--output_dir", type="character", default="./",
	help="output dir [default=%(default)s]", metavar="directory")
parser$add_argument("-p", "--plate", type="character", required=TRUE, action="append",
	help="plate to compare (use multiple times for each)", metavar="group")
parser$add_argument("--zfile_basename", type="character", default="Zscores.csv",
	help="zfile_basename [default=%(default)s]", metavar="Zscores file basename")
# store_true means "int=False unless --int passed, then int=True" (store_false is the inverse)
#parser.add_argument('--int', action='store_true',
#	help='convert values to ints to %(prog)s (default: %(default)s)')
parser$add_argument('--counts', action='store_true',
	help='values are counts not zscores to %(prog)s (default: %(default)s)')
parser$add_argument('--ignore_plate', action='store_true',
	help='ignore the plate and do not include it in the formula to %(prog)s (default: %(default)s)')
opt <- parser$parse_args()


groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

plates=opt$plate
print("Comparing these plates")
print(plates)

owd=opt$output_dir
print("Output dir")
print(owd)
dir.create(owd,showWarnings=F)


library(data.table)

# Multi plate Logistic regression for tile presence on case/control status, adjusting for age, sex, and plate
Z = opt$zscore

#Groups to compare, from the manifest file.
# Order here matters, the first will be coded to 1, the second to 0. So choose the event (aka glioma or pemphigus) to be coded to 1.

date=format(Sys.Date(),"%Y%m%d")

output_base = paste0(owd, "/", gsub(" ","_",
	paste("Multiplate_Peptide_Comparison",
		fs::path_ext_remove(basename(opt$zfile_basename)),
		"type",opt$type,
		paste(groups_to_compare, collapse="-"),
		"Z", Z,"sex",opt$sex,sep="-")))

print(output_base)
dir.create(owd, recursive = TRUE)


# Log the parameter choices into a logfile
logname = paste0(output_base,'.log')

cat("Multi plate Logistic regression for tile presence on case/control status, adjusting for age, sex, and plate.",
	file=logname,sep="\n")
cat(paste0("Z-score threshold = ", Z, "\n"), file = logname, append = TRUE, sep = "\n")
cat("Plates used in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(plates))){
	cat(plates[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)

cat("\nGroups compared in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(groups_to_compare))){
	cat(groups_to_compare[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)



# Read in multiple plate Zscores files.
Zfiles = list()
species_ids = list()

for(i in c(1:length(plates))){
	results=read_zfile(paste(plates[i], opt$zfile_basename, sep="/"))
	species_ids[[i]] = results$species
	Zfiles[[i]] = results$zfile
}
rm(results)

manifest = read_multiple_manifests(plates)

uniq_sub = select_subjects(manifest,opt)
cat(paste0("\nUnique subjects: ", paste(uniq_sub, collapse=",")), file = logname, append = TRUE, sep = "\n")

cat(paste0("\nTotal number of included subjects: ", length(uniq_sub)), file = logname, append = TRUE, sep = "\n")

# Get the overlapping peptides included in each file
common_peps = Reduce(intersect, sapply(species_ids, `[`,1))
print("length(common_peps)")
print(length(common_peps))
print("common_peps[1:5]")
print(common_peps[1:5])
#	[1] "1"     "10"    "100"   "1000"  "10000"

cat(paste0("\nTotal number of included peptides: ", length(common_peps)), file = logname, append = TRUE, sep = "\n")


# Go back into the Z score files and reorder them/cull them to have only the common_peps, in that specific order.
print("for(i in c(1:length(plates))){")
print("IMPORTANT: The column is 'id' NOT 'ids'")
#	[1] "IMPORTANT: The column is 'id' NOT 'ids'"
for(i in c(1:length(plates))){
	print(i)
	print(Zfiles[[i]][1:5,1:5])
	Zfiles[[i]] = Zfiles[[i]][,c("id", common_peps)]
}
#	[1] 1
#	           id   1  10                  100 1000
#	1          id   1  10                  100 1000
#	3    14078-01 0.0 0.0   0.8133400627769187  0.0
#	4 14078-01dup 0.0 0.0 -0.48156984062927694  0.0
#	5    14118-01 0.0 0.0   1.1354027662377901  0.0
#	6 14118-01dup 0.0 0.0    5.141737192497063  0.0



# Convert every Z score pair to a binary 0/1 call based on the Z thresh.


cat("\nStart converting Z scores to peptide binary calls", file = logname, append = TRUE, sep = "\n")


peptide_calls = data.frame(mat.or.vec(length(uniq_sub), length(common_peps)+1))
colnames(peptide_calls) = c("ID", common_peps)
peptide_calls$ID = uniq_sub
print("peptide_calls[1:5,1:5]")
print(peptide_calls[1:9,1:9])
#        ID 1 10 100 1000
#1 14078-01 0  0   0    0
#2 14118-01 0  0   0    0
#3 14127-01 0  0   0    0
#4 14142-01 0  0   0    0
#5 14206-01 0  0   0    0
print(dim(peptide_calls))
#[1]  126 1001


# Loop over every person, and convert all of their peptide Z scores into 0's and 1's.
print("Loop over every person, extract mins and, if z scores, convert all of their peptide Z scores into 0's and 1's.")
print("for(i in c(1:nrow(peptide_calls))){")
for(i in c(1:nrow(peptide_calls))){
	print(i)

	id = peptide_calls$ID[i]
	# Get plate to pull from
	mp = manifest$plate[which(manifest$subject==id)][[1]]

	# Extract the rows containing the duplicate samples, assumes they both contain the id, and that id doesn't match another subjects in a grep (i.e IDs 100 and 1000.)
	rlocs = grep(id, Zfiles[[mp]][,1])

	myZs = data.frame(t(Zfiles[[mp]][c(1,rlocs), which(Zfiles[[mp]][1,] %in% common_peps)]))
	myZs[,c(2:3)] = sapply(myZs[,c(2:3)], as.numeric)

	if(length(which((common_peps == myZs[,1])==FALSE))>0){
		cat(paste("\nError:", id, ":", "plate", mp, ":Peptides out of order, need to ensure they are ordered the same as the common_peps vector. This should never appear.", sep = " "), file=logname, append = TRUE, sep= "\n")
	}
	mymins = apply(myZs[,c(2:3)],1, min)
	if(length(is.na(mymins)) >0){
		mymins[which(is.na(mymins))] = 0
	}

	if( !opt$counts ) {
		peptide_calls[i, -1] = ifelse(mymins > Z, 1, 0)
	} else {
		peptide_calls[i, -1] = mymins
	}

}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(Zfiles)

print("peptide_calls[1:9,1:9]")
print(peptide_calls[1:9,1:9])

cat("peptide_calls[1:20,1:12]", file = logname, append = TRUE, sep = "\n")
cat(capture.output(print(peptide_calls[1:20,1:12])), file = logname, append = TRUE, sep = "\n")

#	Create a shell file for analysis, Leaves the peptide column blank, we will repopulate this with every peptide.
#	I still don't know what a "shell file" is in this context

datfile = build_datfile(uniq_sub,opt)
datfile$peptide = NA

#datfile = data.frame(mat.or.vec(length(uniq_sub),6))
#colnames(datfile) = c("ID", "case", "peptide", "sex", "age", "plate")
#datfile$ID = uniq_sub
#datfile$peptide = NA
#print("for(i in c(1:nrow(datfile))){")
#print("Prep datfile")
#for(i in c(1:nrow(datfile))){
#	print(i)
#	man_loc = which(manifest$subject== datfile$ID[i])[1]
#	datfile$case[i] = ifelse(manifest$group[man_loc] == groups_to_compare[1], 1, 0)
#	datfile$age[i] = manifest$age[man_loc]
#	datfile$sex[i] = manifest$sex[man_loc]
#	datfile$plate[i] = manifest$plate[man_loc]
#}
#datfile$age = as.numeric(datfile$age)
#datfile$sex = as.factor(datfile$sex)
#datfile$plate = as.factor(datfile$plate)
#print(datfile$plate)

print(head(datfile))
cat(capture.output(print(head(datfile))), file = logname, append = TRUE, sep = "\n")

formula = "case ~ peptide"
if( length(unique(datfile$age)) > 1 )
	formula = paste(formula, "age", sep = " + ")
if( length(unique(datfile$sex)) > 1 )
	formula = paste(formula, "sex", sep = " + ")
if( ( length(unique(datfile$plate)) > 1 ) && ( !opt$ignore_plate ) )
	formula = paste(formula, "plate", sep = " + ")

print(paste("formula",formula))
cat(paste("formula",formula), file = logname, append = TRUE, sep = "\n")



# Result File
pvalues = data.frame(mat.or.vec(length(common_peps), 7))
colnames(pvalues) = c("peptide", "species", "freq_case", "freq_control", "beta", "se", "pval")
pvalues$peptide = common_peps

print("pvalues[1:5,1:5]")
print(pvalues[1:5,1:5])

n_case = length(which(datfile$case==1))
n_control = length(which(datfile$case==0))

cat(paste0("\nTotal number of ", groups_to_compare[1], ": ", n_case), file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal number of ", groups_to_compare[2], ": ", n_control), file = logname, append = TRUE, sep = "\n")


# Loop over peptides, populate the datfile, and run analyses.
cat("\nStart loop over peptide logistic regression analysis:", file = logname, append = TRUE, sep = "\n")

print("length(common_peps)")
print(length(common_peps))

cat("length(common_peps)", file = logname, append = TRUE, sep = "\n")
cat(length(common_peps), file = logname, append = TRUE, sep = "\n")

print("for(i in c(1:length(common_peps))){")
for(i in c(1:length(common_peps))){
	print(i)

	# Pull the species assignment of the peptide (this could be more efficient but not really a heavy lookup)
	pvalues$species[i] = species_ids[[1]]$species[which(species_ids[[1]]$id==common_peps[i])][1]

	# Extract all peptide information into the new datfile.
	datfile$peptide = peptide_calls[, which(colnames(peptide_calls)== common_peps[i])]

	if( !opt$counts ) {

		# Calculate the frequency of the peptide in each group
		n_case_pos = length(which(datfile$peptide==1 & datfile$case==1))
		n_control_pos = length(which(datfile$peptide==1 & datfile$case==0))

		pvalues$freq_case[i] = n_case_pos/n_case
		pvalues$freq_control[i] = n_control_pos/n_control

		#	Ignore if all or none of the samples include the peptide. (This is about 70%)
		if( (n_case_pos +n_control_pos) %in% c(0, n_case+n_control)){
			pvalues$beta[i]= NA
			pvalues$se[i] = NA
			pvalues$pval[i] = NA
		}else{
			results = log_reg(datfile,formula,'peptide')
			pvalues$beta[i]= results[1]
			pvalues$se[i] = results[2]
			pvalues$pval[i] = results[3]
		}

	} else {

		#	prep for log transform of 0, or perhaps even negative
		datfile$peptide <- ifelse(datfile$peptide <= 0, 0.001, datfile$peptide)

		#	log transform data to normalize
		datfile$peptide <- log(datfile$peptide)

		pvalues$freq_case[i] = "UNK"
		pvalues$freq_control[i] = "UNK"
		results = log_reg(datfile,formula,'peptide')
		pvalues$beta[i]= results[1]
		pvalues$se[i] = results[2]
		pvalues$pval[i] = results[3]

	}

}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")


print("pvalues[1:5,1:5]")
print(pvalues[1:5,1:5])
#  peptide               species freq_case freq_control         beta
#1       1 Papiine herpesvirus 2       UNK          UNK   0.00240489
#2      10        Vaccinia virus       UNK          UNK   0.31687319
#3     100   Human herpesvirus 3       UNK          UNK  -0.22301611
#4    1000     Hepatitis B virus       UNK          UNK -52.93485789
#5   10000   Human herpesvirus 8       UNK          UNK  -0.69588877


colnames(pvalues) = c("peptide", "species",
	paste0("freq_", groups_to_compare[1]),
	paste0("freq_", groups_to_compare[2]),
	"beta", "se", "pval")

write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],paste0(output_base,'.csv'),
	col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

cat("\n Analysis complete.", file = logname, append = TRUE, sep = "\n")



