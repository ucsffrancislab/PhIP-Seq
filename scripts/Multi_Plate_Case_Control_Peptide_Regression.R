#!/usr/bin/env Rscript

#	Uses a logistic regression framework to assess the association between the presence of each peptide and user specified case control status.
#	Uses a user specifed Z score threshold (i.e. Z >3.5) to call positivity
#	Takes as input a list of directories (e.g. plates), to use in the analysis, the directories should contain the Zscores.t.csv and the manifest* files.
#	Only peptides found in both files are included in the analysis
#	The logistic regression adjusts for age (continuous), sex (M/F factor), and plate (1..n plates, factor)
#	The output (to a user specified directory owd , is a list of peptides, ordered by ascending p-value. It includes peptide ID, species, frequency of peptide in the first group, frequency of the peptide in the second group, and the beta, standard error, and p-value from the multivariable logistic regression (NA if there is no variation in the presence of the peptide overall, e.g. all 0 or all 1).
#	The file will also output a logfile with similar naming convention, which includes the details of the plates used in the analysis, sample sizes, etc.



library("optparse")

option_list = list(
	make_option(c("-a", "--group1"), type="character", default=NULL,
		help="First group to compare", metavar="character"),
	make_option(c("-b", "--group2"), type="character", default=NULL,
		help="Second group to compare", metavar="character"),
	make_option(c("-p", "--plates_to_compare"), type="character", default=NULL,
		help="Comma separated list of plate dirs to compare", metavar="character"),
#	make_option(c("-m", "--manifest"), type="character", default=NULL,
#		help="manifest file name", metavar="character"),
	make_option(c("-o", "--output_dir"), type="character", default="./",
		help="output dir [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#if (is.null(opt$manifest)){
#	print_help(opt_parser)
#	stop("manifest file required.\n", call.=FALSE)
#}

if (is.null(opt$plates_to_compare)){
	print_help(opt_parser)
	stop("plates_to_compare required.\n", call.=FALSE)
}

if (is.null(opt$output_dir)){
	print_help(opt_parser)
	stop("output_dir required.\n", call.=FALSE)
}

if (is.null(opt$group1)){
	print_help(opt_parser)
	stop("group1 required.\n", call.=FALSE)
}

if (is.null(opt$group2)){
	print_help(opt_parser)
	stop("group2 required.\n", call.=FALSE)
}

groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

plates=unlist(strsplit(opt$plates_to_compare, split = ","))
print("Comparing these plates")
print(plates)

owd=opt$output_dir
print("Output dir")
print(owd)
dir.create(owd,showWarnings=F)




# Multi plate Logistic regression for tile presence on case/control status, adjusting for age, sex, and plate
Z= 3.5

# #GBM

# list of paths to the Z score files, each path represents a plate.
#plates= c("/Users/gguerra/Library/CloudStorage/Box-Box/Francis\ _Lab_Share/20241224-Illumina-PhIP/20241224c-PhIP/out.gbm", "/Users/gguerra/Library/CloudStorage/Box-Box/Francis _Lab_Share/20241204-Illumina-PhIP/20241204c-PhIP/out.gbm.test6")
# Directory to pipe all results to
#owd = "/Users/gguerra/Library/CloudStorage/Box-Box/Francis\ _Lab_Share/20250102-PhIP-gbm"
#Groups to compare, from the manifest file.
# Order here matters, the first will be coded to 1, the second to 0. So choose the event (aka glioma or pemphigus) to be coded to 1.
# groups_to_compare = c("case", "control")

#Pemphigus
# plates = c("/Users/gguerra/Library/CloudStorage/Box-Box/Francis _Lab_Share/20241204-Illumina-PhIP/20241204c-PhIP/out.menpem.test6", "/Users/gguerra/Library/CloudStorage/Box-Box/Francis _Lab_Share/20241224-Illumina-PhIP/20241224c-PhIP/out.menpem")
# owd = "/Users/gguerra/Library/CloudStorage/Box-Box/Francis\ _Lab_Share/20250102-PhIP-pemphigus"
# groups_to_compare=c("PF Patient", "Endemic Control" )
# groups_to_compare=c("PF Patient", "Non Endemic Control" )
# groups_to_compare=c("Endemic Control", "Non Endemic Control" )






#date="20250102"
date=format(Sys.Date(),"%Y%m%d")

output_base = paste0(owd, "/", gsub(" ","_",
	paste(date, "Multiplate_Peptide_Comparison", paste(groups_to_compare, collapse="-"),"Prop_test_results", Z,sep="-")))

# Log the parameter choices into a logfile
#logname = paste0(owd, "/", date, "_Multiplate_Peptide_Comparison_", groups_to_compare[1], "_", groups_to_compare[2],"_Prop_test_results_", Z, ".log")
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



# Read in multiple plate Zscores.t files.
Zfiles = list()
species_ids = list()

for(i in c(1:length(plates))){
	Zfilename = paste0(plates[i], "/Zscores.t.csv")
	Zfile= read.csv(Zfilename, sep = ",", header=FALSE)
	Zfile = data.frame(t(Zfile))

	# If in the format of subject, type species, remove subject and type, and remove second row.
	if("subject" %in% Zfile[2,c(1:3)]){
		to_remove= which(Zfile[2,c(1:3)]== "subject")
		Zfile = Zfile[,-to_remove]
	}
	if("type" %in% Zfile[2,c(1:3)]){
		to_remove= which(Zfile[2,c(1:3)]== "type")
		Zfile = Zfile[,-to_remove]
	}

	# Extract the peptide information
	species_id = data.frame(t(Zfile[c(1:2),]))
	colnames(species_id) = species_id[1,]
	species_id = species_id[-1,]
	species_ids[[i]] = species_id
	Zfile = Zfile[-2,]
	colnames(Zfile) = Zfile[1,]
	Zfiles[[i]] = Zfile
	colnames(species_ids[[i]]) = c("id", "species")
}
rm(Zfile)
rm(species_id)


# Read in the multiple manifest files.
mfs  = list()

for(i in c(1:length(plates))){
	# Find the manifest file for the given plate. Requires only ONE manifest file per plate folder
	mfname = list.files(plates[i], pattern="manifest", full.names=TRUE)
	if(length(mfname)!=1){
		print(paste0(plates[i], " needs a single manifest file!"))
	}

	# read in the manifest file
	#mf = read.csv(paste(mfname, sep = ""), sep= ",", header = TRUE)
	mf = read.csv(mfname, sep= ",", header = TRUE)
	# Create a categorical variable, assign all of these the same number to indicate plate.
	mf$plate = i
	mfs[[i]] = mf
}

# Create an aggregate metadata file. This requires identical column structure in the files.

manifest = Reduce(rbind, mfs)
# can get rid of mfs list.
rm(mfs)

# Identify the unique subjects to include in the analyses.

uniq_sub = unique(manifest$subject[which(manifest$group %in% groups_to_compare)])

cat(paste0("\nTotal number of included subjects: ", length(uniq_sub)), file = logname, append = TRUE, sep = "\n")

# Get the overlapping peptides included in each file
common_peps = Reduce(intersect, sapply(species_ids, `[`,1))

cat(paste0("\nTotal number of included peptides: ", length(common_peps)), file = logname, append = TRUE, sep = "\n")


# Go back into the Z score files and reorder them/cull them to have only the common_peps, in that specific order.
for(i in c(1:length(plates))){
	Zfiles[[i]] = Zfiles[[i]][,c("id", common_peps)]
}

# Convert every Z score pair to a binary 0/1 call based on the Z thresh.

cat("\nStart converting Z scores to peptide binary calls", file = logname, append = TRUE, sep = "\n")

peptide_calls = data.frame(mat.or.vec(length(uniq_sub), length(common_peps)+1))
colnames(peptide_calls) = c("ID", common_peps)
peptide_calls$ID = uniq_sub

# Loop over every person, and convert all of their peptide Z scores into 0's and 1's.
for(i in c(1:nrow(peptide_calls))){

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

	peptide_calls[i, -1] = ifelse(mymins > Z, 1, 0)
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(Zfiles)

# Create a shell file for analysis, Leaves the peptide column blank, we will repopulate this with every peptide.

datfile = data.frame(mat.or.vec(length(uniq_sub),6))
colnames(datfile) = c("ID", "case", "peptide", "sex", "age", "plate")
datfile$ID = uniq_sub
datfile$peptide = NA
for(i in c(1:nrow(datfile))){
	man_loc = which(manifest$subject== datfile$ID[i])[1]
	datfile$case[i] = ifelse(manifest$group[man_loc] == groups_to_compare[1], 1, 0)
	datfile$age[i] = manifest$age[man_loc]
	datfile$sex[i] = manifest$sex[man_loc]
	datfile$plate[i] = manifest$plate[man_loc]
}


datfile$age = as.numeric(datfile$age)
datfile$sex = as.factor(datfile$sex)
datfile$plate = as.factor(datfile$plate)


#----- Shell function for logistic regression analysis.
log_reg = function(df){


	# A simple model that simply adjusts for plate/batch in the model. When the number of plates becomes large, a mixed effects regression model should be considered.
	# as there are likely differences in the peptide calling sensitivity between plates, and so peptide probably has different associations with case based on plate.
	logitmodel = "case~ peptide +age + sex + plate"

	logit_fun = glm(as.formula(logitmodel), data = df, family=binomial(link="logit"))

	go= summary(logit_fun)
	beta = go$coefficients[2,1]
	se = go$coefficients[2,2]
	pval = go$coefficients[2,4]
	return(c(beta, se, pval))

}
#-------


# Result File
pvalues = data.frame(mat.or.vec(length(common_peps), 7))
colnames(pvalues) = c("peptide", "species",  "freq_case", "freq_control", "beta", "se", "pval")
pvalues$peptide = common_peps

n_case = length(which(datfile$case==1))
n_control = length(which(datfile$case==0))

cat(paste0("\nTotal number of ", groups_to_compare[1], ": ", n_case), file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal number of ", groups_to_compare[2], ": ", n_control), file = logname, append = TRUE, sep = "\n")


# Loop over peptides, populate the datfile, and run analyses.
cat("\nStart loop over peptide logistic regression analysis:", file = logname, append = TRUE, sep = "\n")

for(i in c(1:length(common_peps))){

	# Pull the species assignment of the peptide (this could be more efficient but not really a heavy lookup)
	pvalues$species[i] = species_ids[[1]]$species[which(species_ids[[1]]$id==common_peps[i])][1]

	# Extract all peptide information into the new datfile.
	datfile$peptide = peptide_calls[, which(colnames(peptide_calls)== common_peps[i])]

	# Calculate the frequency of the peptide in each group
	n_case_pos = length(which(datfile$peptide==1 & datfile$case==1))
	n_control_pos = length(which(datfile$peptide==1 & datfile$case==0))

	pvalues$freq_case[i] = n_case_pos/n_case
	pvalues$freq_control[i] = n_control_pos/n_control

	if( (n_case_pos +n_control_pos) %in% c(0, n_case+n_control)){
		pvalues$beta[i]= NA
		pvalues$se[i] = NA
		pvalues$pval[i] = NA
	}else{
		results =log_reg(datfile)
		pvalues$beta[i]= results[1]
		pvalues$se[i] = results[2]
		pvalues$pval[i] = results[3]
	}
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")


colnames(pvalues) = c("peptide", "species",
	paste0("freq_", groups_to_compare[1]),
	paste0("freq_", groups_to_compare[2]),
	"beta", "se", "pval")

write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],paste0(output_base,'.csv'),
	col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

cat("\n Analysis complete.", file = logname, append = TRUE, sep = "\n")



