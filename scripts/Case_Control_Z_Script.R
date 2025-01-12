#!/usr/bin/env Rscript

#	Assesses differences between epitopes for two groups.
#	For a user provided Z-score threshold, and two groups (chosen from the manifest "type" column), this will compute a p-value for the difference in proportion of samples with that tile present between the two groups. Writes a file in the same data directory, with name "Tile Comparison*" (this naming is new so it doesnt match the output names in the older test folders.
#	Output p-values are not adjusted for any multiple testing. NA p values are produced when tile proportion is 1 in both or 0 in both.





#library("optparse")
#
#option_list = list(
#  make_option(c("-z", "--zscore"), type="double", default=3.5,
#    help="Zscore threshold", metavar="character"),
#	make_option(c("-a", "--group1"), type="character", default=NULL,
#		help="First group to compare", metavar="character"),
#	make_option(c("-b", "--group2"), type="character", default=NULL,
#		help="Second group to compare", metavar="character"),
#	make_option(c("-g", "--groups_to_compare"), type="character", default=NULL,
#		help="Comma separated list of groups to compare", metavar="character"),
#	make_option(c("-m", "--manifest"), type="character", default=NULL,
#		help="manifest file name", metavar="character"),
#	make_option(c("-d", "--working_dir"), type="character", default="./",
#		help="working dir [default= %default]", metavar="character")
#);
#
#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);
#
#
#if (is.null(opt$manifest)){
#	print_help(opt_parser)
#	stop("manifest file required.\n", call.=FALSE)
#}
#
#if (is.null(opt$group1)){
#	print_help(opt_parser)
#	stop("group1 required.\n", call.=FALSE)
#}
#
#if (is.null(opt$group2)){
#	print_help(opt_parser)
#	stop("group2 required.\n", call.=FALSE)
#}



library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
parser <- ArgumentParser(description=scriptname)
parser$add_argument("-a", "--group1", type="character", required=TRUE,
	help="first group to compare", metavar="group")
parser$add_argument("-b", "--group2", type="character", required=TRUE,
	help="second group to compare", metavar="group")
parser$add_argument("-z", "--zscore", type="double", default=3.5,
	help="zscore threshold", metavar="double")
parser$add_argument("-m", "--manifest", type="character", default=NULL, required=TRUE,
	help="manifest file name", metavar="manifest")
parser$add_argument("-d", "--working_dir", type="character", default="./",
	help="working dir [default=%(default)s]", metavar="directory")
opt <- parser$parse_args()




# For each tile, compares proportion of present in each group using 2-prop test. Reports proportions and associated p-value.

# Input parameters
#groups_to_compare = c("case", "control")
#groups_to_compare=c("PF Patient", "Endemic Control" )
#groups_to_compare=unlist(strsplit(opt$groups_to_compare, split = ","))
#groups_to_compare=unlist(strsplit(opt$groups, split = ","))

groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

#Z = 3.5
Z = opt$zscore




# Read in the Z-score file  (wihtout transpose and remove the transpose line)

Zfile = read.csv(paste(opt$working_dir, "Zscores.t.csv", sep = "/"), sep = ",", header=FALSE)

Zfile = data.frame(t(Zfile))

# Read in the metadata file

meta = read.csv(opt$manifest, sep= ",", header = TRUE)


# If in the format of subject, type species, remove subject and type, and remove second row.
if("subject" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "subject")
	Zfile = Zfile[,-to_remove]
}

if("type" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "type")
	Zfile = Zfile[,-to_remove]
}

print("Extract the peptide information")

species_id = data.frame(t(Zfile[c(1:2),]))
colnames(species_id) = species_id[1,]
species_id = species_id[-1,]

Zfile = Zfile[-2,]

print("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])
to_keep = 1
for(u in uniqid){
	possible_ids = grep(u, Zfile[,1])
	mids = Zfile[possible_ids,1]
	locs = grepl("dup", mids)
	to_keep = c(to_keep,possible_ids[c(which(locs==FALSE)[1],which(locs==TRUE)[1]) ] )
}

Zfile1 = Zfile[to_keep,]
rm(Zfile)

print("Create a shell file for analysis")

datfile = data.frame(mat.or.vec(length(unique(Zfile1[,1]))-1,3))
colnames(datfile) = c("ID", "case", "peptide")
datfile$ID = uniqid
for(i in c(1:nrow(datfile))){
	datfile$case[i] = meta$group[which(meta$subject== datfile$ID[i])[1]]
}


# Result File
pvalues = data.frame(mat.or.vec(ncol(Zfile1)-1, 5))
colnames(pvalues) = c("peptide", "species", "freq_case", "freq_control", "pval")


print("Pick a peptide:")
pep_index = 1
for(pep_index in c(1:(ncol(Zfile1)-1))){
	pepcol = pep_index+1
	peptide = Zfile1[1,(pepcol)]
	pvalues$peptide[pep_index] = peptide

	pvalues$species[pep_index] = species_id[which(species_id$id==peptide), 2]

	print(peptide)

	datfile$peptide = NA
	for(i in c(1:nrow(datfile))){
		Zs = Zfile1[which(Zfile1[,1] == datfile$ID[i]), pepcol]
		datfile$peptide[i] = ifelse(min(as.numeric(Zs) >Z), 1,0)
	}

	dfl = datfile[complete.cases(datfile), ]
	pvalues$pval[pep_index] = NA
	if(length(which(dfl$case==groups_to_compare[1]))>4 & length(which(dfl$case==groups_to_compare[2]))>4){
		# If its absent in everyone, report that
		if(sum(dfl$peptide)==0 ){
			pvalues$freq_case[pep_index] = 0
			pvalues$freq_control[pep_index] = 0
			pvalues$pval[pep_index] = NA
		}else{

			# If its present in everyone, report that
			if(sum(dfl$peptide)==nrow(dfl) ){
				pvalues$freq_case[pep_index] = 1
				pvalues$freq_control[pep_index] = 1
				pvalues$pval[pep_index] = NA
			}else{

				# If its not present in everyone, do some analysis
				# For now, a proportion test
				n_cases = length(which(dfl$case==groups_to_compare[1]))
				n_control = length(which(dfl$case==groups_to_compare[2]))
				n_case_success = length(which(dfl$case==groups_to_compare[1] & dfl$peptide ==1))
				n_control_success = length(which(dfl$case==groups_to_compare[2] & dfl$peptide==1))

				prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control), p = NULL, alternative = "two.sided", correct = TRUE)
				pvalues$freq_case[pep_index] = n_case_success/n_cases
				pvalues$freq_control[pep_index] = n_control_success/n_control
				pvalues$pval[pep_index] = prop$p.value

			}
		}
	}
}

print("close loop over peptides")

colnames(pvalues) = c("peptide", "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")


outfile=paste0(opt$working_dir, "/",
	gsub(" ","_",paste("Tile_Comparison", paste(groups_to_compare[1:2],collapse="-"), "Prop_test_results", Z, sep="-")), ".csv")


print(paste0("Writing ",outfile))
write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),], outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


