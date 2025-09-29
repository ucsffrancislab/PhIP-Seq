#!/usr/bin/env Rscript

#	Uses a logistic regression framework to assess the association between the presence of each virus and user specified case control status.
#	Uses roughly VirScan's definition of seropositivity (Presence of a public epitope and >1 unique reads hit)
#	Takes as input a list of directories (e.g. plates), to use in the analysis, the directories should contain the seropositive.csv and the manifest* files.
#	Only viruses found in both files are included in the analysis
#	Only uses the "Before" datatype, where the public epitopes are checked prior to filtering for unique reads.
#	The logistic regression adjusts for age (continuous), sex (M/F factor), and plate (1..n plates, factor)
#	The output (to a user specified directory owd , is a list of viral species, ordered by ascending p-value. It includes species, frequency of virus in the first group, frequency of the virus in the second group, and the beta, standard error, and p-value from the multivariable logistic regression (NA if there is no variation in the presence of the virus overall, e.g. all 0 or all 1).
#	The file will also output a logfile with similar naming convention, which includes the details of the plates used in the analysis, sample sizes, etc.

options(warn = 1)

library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
scriptdir=dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(paste(scriptdir,'GenoLib.R',sep='/'))
parser <- ArgumentParser(description=scriptname)
parser$add_argument("--study", type="character", default="",
	help="limit study", metavar="study")
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
parser$add_argument("-p", "--plate", type="character", required=TRUE, action="append",
	help="plate to compare (use multiple times for each)", metavar="group")
parser$add_argument("-o", "--output_dir", type="character", default="./",
	help="output dir [default=%(default)s]", metavar="directory")
parser$add_argument("--sfile_basename", type="character", default="seropositive.csv",
	help="sfile_basename [default=%(default)s]", metavar="seropositive file basename")
#parser$add_argument("--keep_only_B", action="store_true",
#	help="Keep only those ids with '_B'")
parser$add_argument("--keep_all_ids", action="store_true",
	help="Keep all ids. Not just '_B'")
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

# Multi plate Logistic regression for seropositivity (VirScan calls) on case/control status, adjusting for age, sex, and plate

# # list of paths to the Z score files, each path represents a plate.
# #Groups to compare, from the manifest file.
# # Order here matters, the first will be coded to 1, the second to 0. So choose the event (aka glioma or pemphigus) to be coded to 1.
# groups_to_compare = c("case", "control")

date=format(Sys.Date(),"%Y%m%d")

Z = opt$zscore

#	paste(date, "Multiplate_VirScan_Seropositivity_Comparison",
output_base = paste0(owd, "/", gsub(" ","_",
	paste("Multiplate_VirScan_Seropositivity_Comparison",
	fs::path_ext_remove(basename(opt$sfile_basename)),
	opt$study, opt$type, paste(groups_to_compare, collapse="-"),
	"Z",Z,
	"sex",opt$sex, sep="-")))

# Log the parameter choices into a logfile
logname = paste0(output_base,'.log')

cat("Multi plate Logistic regression for presence of virus on case/control status, adjusting for age, sex, and plate. Using VirScan's parameters for virus calling.",file=logname,sep="\n")

cat("\nPlates used in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(plates))){
	cat(plates[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)

cat("\nGroups compared in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(groups_to_compare))){
	cat(groups_to_compare[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)

posfiles = list()

# Read in multiple plate seropositivity files.
for(i in c(1:length(plates))){

	posfile <- data.frame(data.table::fread(paste(plates[i],opt$sfile_basename,sep="/"),
		sep = ",", header=TRUE), check.names = FALSE )

	if( opt$keep_all_ids ){
		posfile1 = posfile
	}else{
		posfile1 = posfile[grep("_B$", posfile$id), ]
	}
	rm(posfile)

	posfiles[[i]] = posfile1

}# close loop over plates.

manifest = read_multiple_manifests(plates)

uniq_sub = select_subjects(manifest,opt)
cat(paste0("\nTotal number of included subjects: ", length(uniq_sub)), file = logname, append = TRUE, sep = "\n")

# Identify the viruses that are in both files, and subset each posfile to that, with maintained order.
common_virs = Reduce(intersect, lapply(posfiles,function(x) colnames(x)))

for(i in c(1:length(plates))){
	posfiles[[i]] = posfiles[[i]][,common_virs]
}
common_virs = common_virs[-c(1:3)]

cat(paste0("\nTotal number of included viruses: ", length(common_virs)), file = logname, append = TRUE, sep = "\n")

# Convert every virus call into a binary 0/1 call based on >1 hits (Jake already filtered for public epitopes).

cat("\nStart converting viral calls to binary calls", file = logname, append = TRUE, sep = "\n")

viral_calls = data.frame(mat.or.vec(length(uniq_sub), length(common_virs)+1))
colnames(viral_calls) = c("ID", common_virs)
viral_calls$ID = uniq_sub

# Loop over every person, and convert all of their peptide Z scores into 0's and 1's.
for(i in c(1:nrow(viral_calls))){
	print(paste("Looping:",i,":",nrow(viral_calls)))

	id = viral_calls$ID[i]
	# Get plate to pull from
	mp = manifest$plate[which(manifest$subject==id)][[1]]

	# Extract the rows containing the duplicate samples, assumes they both contain the id,
	#	and that id doesn't match another subjects in a grep (i.e IDs 100 and 1000.)
	rloc = which(posfiles[[mp]][,1] == id)[1]

	myCounts = data.frame(t(posfiles[[mp]][rloc, which(colnames(posfiles[[mp]]) %in% common_virs)]))
	myCounts[,1]=as.numeric(myCounts[,1])
	if(length(which((common_virs == row.names(myCounts))==FALSE))>0){
		cat(paste("\nError:", id, ":", "plate", mp,
			":Viruses out of order, need to ensure they are ordered the same as the common_virs vector. This should never appear.",
			sep = " "), file=logname, append = TRUE, sep= "\n")
	}

	viral_calls[i, -1] = ifelse(myCounts > 1, 1, 0)
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(posfiles)


# Create a shell file for analysis, Leaves the virus column blank, we will repopulate this with every peptide.

datfile = build_datfile(uniq_sub,opt)
datfile$virus = NA

print(head(datfile))
cat(capture.output(print(head(datfile))), file = logname, append = TRUE, sep = "\n")

formula = "case ~ virus"
if( length(unique(datfile$age)) > 1 )
	formula = paste(formula, "age", sep = " + ")
if( length(unique(datfile$sex)) > 1 )
	formula = paste(formula, "sex", sep = " + ")
if( ( length(unique(datfile$plate)) > 1 ) && ( !opt$ignore_plate ) )
	formula = paste(formula, "plate", sep = " + ")

print(paste("formula",formula))
cat(paste("formula",formula), file = logname, append = TRUE, sep = "\n")


# Result File
pvalues = data.frame(mat.or.vec(length(common_virs), 6))
colnames(pvalues) = c( "species", "freq_case", "freq_control", "beta", "se", "pval")
pvalues$species = common_virs

n_case = length(which(datfile$case==1))
n_control = length(which(datfile$case==0))

cat(paste0("\nTotal number of ", groups_to_compare[1], ": ", n_case), file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal number of ", groups_to_compare[2], ": ", n_control), file = logname, append = TRUE, sep = "\n")


# Loop over viruses, populate the datfile, and run analyses.
cat("\nStart loop over virus logistic regression analysis:", file = logname, append = TRUE, sep = "\n")

for(i in c(1:length(common_virs))){
	print(paste("Looping:",i,":",length(common_virs)))

	# Extract all virus information into the new datfile.
	datfile$virus = viral_calls[, which(colnames(viral_calls)== common_virs[i])]

	# Calculate the frequency of the peptide in each group
	n_case_pos = length(which(datfile$virus==1 & datfile$case==1))
	n_control_pos = length(which(datfile$virus==1 & datfile$case==0))

	pvalues$freq_case[i] = n_case_pos/n_case
	pvalues$freq_control[i] = n_control_pos/n_control

	if( (n_case_pos +n_control_pos) %in% c(0, n_case+n_control)){
		pvalues$beta[i]= NA
		pvalues$se[i] = NA
		pvalues$pval[i] = NA
	}else{
		results = log_reg(datfile,formula,'virus')
		pvalues$beta[i]= results[1]
		pvalues$se[i] = results[2]
		pvalues$pval[i] = results[3]
	}

}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

colnames(pvalues) = c( "species",
	paste0("freq_", groups_to_compare[1]),
	paste0("freq_", groups_to_compare[2]),
	"beta", "se", "pval")

write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],paste0(output_base,'.csv'),
	col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE )

cat("\n Analysis complete.", file = logname, append = TRUE, sep = "\n")


