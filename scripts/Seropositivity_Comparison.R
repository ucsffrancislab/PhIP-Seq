#!/usr/bin/env Rscript

#	Uses the VirScan defined viral calls to assess difference in seropositivity between two user defined groups.
#	Specifically, it only utilizes the "seropositive.csv" file to make calls.
#	Currently it is hardcoded to call a virus positive the value in the seropositive.csv file is >1.
#		This may not be exact VirScan threshold, but I recall their thresholds were basically this.
#	Only uses the "_B" columns, where the Public epitopes were assessed before viral scoring.
#	Only makes calls for viruses with known public epitopes.
#	Outputs a file "Seropositivity_Prop*" which also indicates the two groups compared.

options(warn = 1)

library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
scriptdir=dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(paste(scriptdir,'GenoLib.R',sep='/'))
parser <- ArgumentParser(description=scriptname)
parser$add_argument("-s", "--sex", type="character", default="",
	help="limit sex", metavar="sex")
parser$add_argument("-t", "--type", type="character", default="",
	help="limit type", metavar="type")
parser$add_argument("-a", "--group1", type="character", required=TRUE,
	help="first group to compare", metavar="group")
parser$add_argument("-b", "--group2", type="character", required=TRUE,
	help="second group to compare", metavar="group")
parser$add_argument("-z", "--zscore", type="double", default=3.5,
	help="zscore threshold", metavar="double")
parser$add_argument("-m", "--manifest", type="character", default=NULL, required=TRUE,
	help="manifest file name", metavar="manifest")
parser$add_argument("--output_dir", type="character", default="./",
	help="output dir [default=%(default)s]", metavar="directory")
parser$add_argument("--sfilename", type="character", default="out/seropositive.csv",
	help="seropositive filename [default=%(default)s]", metavar="Zscores file")
#parser$add_argument("-d", "--working_dir", type="character", default="./",
#	help="working dir [default=%(default)s]", metavar="directory")
#parser$add_argument("--keep_only_B", action="store_true",
#	help="Keep only those ids with '_B'")
parser$add_argument("--keep_all_ids", action="store_true",
	help="Keep all ids. Not just '_B'")
opt <- parser$parse_args()


Z = opt$zscore

library(data.table)

# Seropositivity Test

# Input parameters
groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

posfile <- data.frame(data.table::fread(opt$sfilename, sep = ",", header=TRUE))

print("Keep only 'Before' samples")

if( opt$keep_all_ids ){
	posfile1 = posfile
}else{
	posfile1 = posfile[grep("_B$", posfile$id), ]
}
rm(posfile)

print("Read in the manifest file")
manifest <- data.frame(data.table::fread(opt$manifest, sep = ",", header=TRUE))

print("Unique samples to keep")
uniq_sub = select_subjects(manifest,opt)

print("Keep only the first occurrence of each sample name")
posfile2 =  posfile1[which(posfile1$subject %in% uniq_sub),]
rm(posfile1)

pvalues = data.frame(mat.or.vec(ncol(posfile2)-3, 4))
colnames(pvalues) = c("species", "freq_case", "freq_control", "pval")

cases = unique(manifest$subject[which(manifest$group %in% groups_to_compare[1])])
cases = intersect(cases,uniq_sub)
print(paste("length(cases) :",length(cases)))
controls = unique(manifest$subject[which(manifest$group %in% groups_to_compare[2])])
controls = intersect(controls,uniq_sub)
print(paste("length(controls) :",length(controls)))

df_colnames = colnames(posfile2)
for(sp in c(1:(nrow(pvalues)))){
	print(paste("Looping:",sp,":",nrow(pvalues)))
	species = df_colnames[sp+3]
	pvalues$species[sp] = species
	n_cases = length(cases)
	n_control = length(controls)
	n_case_success = length(which(posfile2$subject %in% cases & posfile2[,sp+3] >1))
	n_control_success = length(which(posfile2$subject %in% controls  & posfile2[,sp+3]>1))

	pvalues$freq_case[sp] = n_case_success/n_cases
	pvalues$freq_control[sp] = n_control_success/n_control
	prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control),
		p = NULL, alternative = "two.sided", correct = TRUE)
	pvalues$pval[sp] = prop$p.value
}

colnames(pvalues) = c( "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]


outfile = paste0(opt$output_dir, "/",
	gsub(" ","-",paste("Seropositivity_Prop",
		fs::path_ext_remove(basename(opt$sfilename)),
		opt$type, paste(groups_to_compare[1:2],collapse="-"),
		"Z",Z,"sex",opt$sex,sep="-")), ".csv")


print(paste0("Writing ",outfile))
write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],
	outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

