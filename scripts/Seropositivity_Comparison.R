#!/usr/bin/env Rscript

#	Uses the VirScan defined viral calls to assess difference in seropositivity between two user defined groups.
#	Specifically, it only utilizes the "seropositive.csv" file to make calls.
#	Currently it is hardcoded to call a virus positive the value in the seropositive.csv file is >1. This may not be exact VirScan threshold, but I recall their thresholds were basically this.
#	Only uses the "_B" columns, where the Public epitopes were assessed before viral scoring.
#	Only makes calls for viruses with known public epitopes.
#	Outputs a file "Seropositivity_Prop_test_results*" which also indicates the two groups compared.


library("optparse")

option_list = list(
  make_option(c("-z", "--zscore"), type="double", default=3.5,
    help="Zscore threshold", metavar="character"),
	make_option(c("-a", "--group1"), type="character", default=NULL,
		help="First group to compare", metavar="character"),
	make_option(c("-b", "--group2"), type="character", default=NULL,
		help="Second group to compare", metavar="character"),
#	make_option(c("-g", "--groups_to_compare"), type="character", default=NULL,
#		help="Comma separated list of groups to compare", metavar="character"),
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

if (is.null(opt$group1)){
	print_help(opt_parser)
	stop("group1 required.\n", call.=FALSE)
}

if (is.null(opt$group2)){
	print_help(opt_parser)
	stop("group2 required.\n", call.=FALSE)
}

Z=opt$zscore


# Seropositivity Test

# Input parameters

#groups_to_compare = c("case", "control" )
#groups_to_compare=unlist(strsplit(opt$groups_to_compare, split = ","))

groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

posfile = read.csv(paste0(opt$working_dir, paste("/seropositive",Z,"csv",sep=".")), header = TRUE, sep = ",")

print("Keep only 'Before' samples")

posfile1 = posfile[grep("_B", posfile$id), ]
rm(posfile)

print("Read in the metadata file")

meta = read.csv(opt$manifest, sep= ",", header = TRUE)

print("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])

print("Keep only the first occurrence of each sample name")
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

colnames(pvalues) = c( "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]

outfile=paste0(opt$working_dir, "/",
	gsub(" ","-",paste("Seropositivity_Prop_test_results", paste(groups_to_compare[1:2],collapse="-"), sep="-")), ".csv")


print(paste0("Writing ",outfile))
write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),], outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

