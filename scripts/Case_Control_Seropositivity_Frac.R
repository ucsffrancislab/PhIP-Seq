#!/usr/bin/env Rscript

#	Assesses differences in the proportion of samples with virus called present between two user specified groups (from the manifest "type" column).
#	Uses our own definition of seropositivity (Virus present if >"Vir_frac" proportion of tiles are hit at Z > "Z_thresh"). I typically have used Vir_frac = 0.02 or 0.05 based on eyeballing the Viral_Frac_Hits table generated above.
#	Outputs a file in the same test directory called "Viral_Sero_test_results*" with indicators of the groups and parameters used.



library("optparse")

option_list = list(
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


# For each virus, make a call of positive or negative based on at least 5% of possible tiles hitting, then measure proportion.
# Input parameters

#groups_to_compare = c("case", "control" )
#groups_to_compare=c("PF Patient", "Endemic Control" )
#groups_to_compare=unlist(strsplit(opt$groups_to_compare, split = ","))
#groups_to_compare=unlist(strsplit(opt$groups, split = ","))

groups_to_compare = c(opt$group1,opt$group2)

print("Comparing these groups")
print(groups_to_compare)


Z_thresh = 3.5
Vir_frac = 0.05
virfracfilename = paste0("Viral_Frac_Hits_Z_",Z_thresh,".csv")

print("Read in the metadata file")

meta = read.csv( opt$manifest, sep= ",", header = TRUE)

print("Read in the VirFrac file")
vir_fracs = read.csv(paste(opt$working_dir, virfracfilename, sep = "/"), header = TRUE)

# # Read in the viral score file
# vs = read.csv(paste(mwd, "/", virfilename, sep =""),sep = ",", header= FALSE )
# vir_score = data.frame(t(vs))
# colnames(vir_score) = vir_score[1,]
# vir_score = vir_score[-1,]


print("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])
#vir_score = vir_score[which(vir_score$id %in% uniqid),]
vir_fracs = vir_fracs[which(vir_fracs$id %in% uniqid),]

cases = unique(meta$subject[which(meta$group %in% groups_to_compare[1])])
controls = unique(meta$subject[which(meta$group %in% groups_to_compare[2])])

print("Create a shell file for analysis")
pvalues = data.frame(mat.or.vec(ncol(vir_fracs)-1, 4))
colnames(pvalues) = c( "species", "freq_case", "freq_control", "pval")
pvalues$species = colnames(vir_fracs)[-1]


for(i in c(1:nrow(pvalues))){
	species = pvalues$species[i]
	n_cases = length(cases)
	n_control = length(controls)

	# For each case, determine number of ids that are >3.5 in both reps

	n_case_success = length(which(as.numeric(vir_fracs[which(vir_fracs$id %in% cases), which(colnames(vir_fracs) == species)]) > Vir_frac))
	n_control_success = length(which(as.numeric(vir_fracs[which(vir_fracs$id %in% controls), which(colnames(vir_fracs) == species)]) > Vir_frac))

	prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control), p = NULL, alternative = "two.sided", correct = TRUE)
	pvalues$freq_case[i] = n_case_success/n_cases
	pvalues$freq_control[i] = n_control_success/n_control
	pvalues$pval[i] = prop$p.value
}

colnames(pvalues) = c( "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]


outfile=paste0(opt$working_dir, "/",
	gsub(" ","_", paste("Viral_Sero_test_results", paste(groups_to_compare[1:2],collapse="-"),
		"Vir_hit_frac", Vir_frac, "Z", Z_thresh, sep="-")), ".csv")

print(paste0("Writing ",outfile))
write.table(opvalues, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


