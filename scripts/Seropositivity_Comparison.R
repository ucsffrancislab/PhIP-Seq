#!/usr/bin/env Rscript

#	Uses the VirScan defined viral calls to assess difference in seropositivity between two user defined groups.
#	Specifically, it only utilizes the "seropositive.csv" file to make calls.
#	Currently it is hardcoded to call a virus positive the value in the seropositive.csv file is >1. This may not be exact VirScan threshold, but I recall their thresholds were basically this.
#	Only uses the "_B" columns, where the Public epitopes were assessed before viral scoring.
#	Only makes calls for viruses with known public epitopes.
#	Outputs a file "Seropositivity_Prop_test_results*" which also indicates the two groups compared.


library("optparse")

option_list = list(
	make_option(c("-m", "--manifest"), type="character", default=NULL,
		help="manifest file name", metavar="character"),
	make_option(c("-d", "--working_dir"), type="character", default="./",
		help="working dir [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$manifest)){
	print_help(opt_parser)
	stop("manifest file required.\n", call.=FALSE)
}



# Seropositivity Test

# Input parameters
plate = "gbm"
groups_to_compare = c("case", "control" )

resdir = paste("out.", plate, ".test5", sep = "")

posfile = read.csv(paste(opt$working_dir, "/", resdir, "/seropositive.csv", sep =""), header = TRUE, sep = ",")

message("Keep only 'Before' samples")

posfile1 = posfile[grep("_B", posfile$id), ]
rm(posfile)

message("Read in the metadata file")

meta = read.csv(opt$manifest, sep= ",", header = TRUE)

message("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])

message("Keep only the first occurrence of each sample name")
posfile2 =  posfile1[which(posfile1$subject %in% uniqid),]
rm(posfile1)

pvalues = data.frame(mat.or.vec(ncol(posfile2)-3, 4))
colnames(pvalues) = c("species", "freq_case", "freq_control", "pval")

cases =unique(meta$subject[which(meta$group %in% groups_to_compare[1])])
controls =unique(meta$subject[which(meta$group %in% groups_to_compare[2])])

df_colnames = colnames(posfile2)
for(sp in c(1:(nrow(pvalues)))){
	species = df_colnames[sp+3]
	pvalues$species[sp] = species
	n_cases = length(cases)
	n_control = length(controls)
	n_case_success = length(which(posfile2$subject %in% cases & posfile2[,sp+3] >1))
	n_control_success = length(which(posfile2$subject %in% controls  & posfile2[,sp+3]>1))

	pvalues$freq_case[sp] = n_case_success/n_cases
	pvalues$freq_control[sp] = n_control_success/n_control
	prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control), p = NULL, alternative = "two.sided", correct = TRUE)
	pvalues$pval[sp] = prop$p.value
}

colnames(pvalues) = c( "species", paste("freq_", groups_to_compare[1], sep= ""), paste("freq_", groups_to_compare[2], sep= ""), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]

write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),], paste(opt$working_dir, "/", resdir, "/Seropositivity_Prop_test_results_", groups_to_compare[1], "_", groups_to_compare[2], ".csv", sep = ""), col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

