#!/usr/bin/env Rscript

#	Assesses differences between epitopes for two groups.
#	For a user provided Z-score threshold, and two groups (chosen from the manifest "type" column),
#	this will compute a p-value for the difference in proportion of samples with that tile present
#	between the two groups. Writes a file in the same data directory, with name "Tile Comparison*"
#		(this naming is new so it doesnt match the output names in the older test folders.
#	Output p-values are not adjusted for any multiple testing.
#		NA p values are produced when tile proportion is 1 in both or 0 in both.


library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
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
parser$add_argument("--zfilename", type="character", default="out/Zscores.csv",
	help="zfilename [default=%(default)s]", metavar="Zscores file")
#parser$add_argument("-d", "--working_dir", type="character", default="./",
#	help="working dir [default=%(default)s]", metavar="directory")
opt <- parser$parse_args()


# For each tile, compares proportion of present in each group using 2-prop test. Reports proportions and associated p-value.

# Input parameters

groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

Z = opt$zscore

library(data.table)

# Read in the Z-score file  (without transpose and remove the transpose line)
Zfile <- data.frame(data.table::fread(opt$zfilename, sep = ",", header=FALSE))
print(Zfile[1:5,1:5])

# Read in the metadata file
meta <- data.frame(data.table::fread(opt$manifest, sep = ",", header=TRUE))


# If in the format of subject, type species, remove subject and type, and remove second row.
if("subject" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "subject")
	Zfile = Zfile[,-to_remove]
}
#	Really should just keep the subject and drop the id column.


if("type" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "type")
	Zfile = Zfile[,-to_remove]
}
print(Zfile[1:5,1:5])

print("Extract the peptide information")

species_id = data.frame(t(Zfile[c(1:2),]))
colnames(species_id) = species_id[1,]
species_id = species_id[-1,]

Zfile = Zfile[-2,]

print("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])
print("length(uniqid)")
print(length(uniqid))

print("uniqid")
print(uniqid)

print("dim(Zfile)")
print(dim(Zfile))


#	something is off here
#	to_keep is a list of column index ids? 
#	grepping for an id has the side effect of getting more than you actually want
#	[1] "4439"
#	[1] 26 27 74 75
#	[1] "1443901"    "1443901dup" "4439"       "4439dup"   
#	[1] FALSE  TRUE FALSE  TRUE

to_keep = 1
for(u in uniqid){

	possible_ids = grep(paste0("^",u,"$|^",u,"dup$"), Zfile[,1])

	mids = Zfile[possible_ids,1]
	locs = grepl("dup", mids)
	to_keep = c(to_keep,possible_ids[c(which(locs==FALSE)[1],which(locs==TRUE)[1]) ] )
}



print("to_keep")
print(to_keep)

print("length(to_keep)")
print(length(to_keep))

Zfile1 = Zfile[to_keep,]
rm(Zfile)

print("dim(Zfile1)")
print(dim(Zfile1))

print("Create a shell file for analysis")

print("Zfile1[,1]")
print(Zfile1[,1])

print("unique(Zfile1[,1])")
print(unique(Zfile1[,1]))

print("length(unique(Zfile1[,1]))")
print(length(unique(Zfile1[,1])))

#	datfile, IMO, is being created bigger than necessary. Because it is a multiple of the length of ids, it works, FOR NOW.
#	gonna leave it alone

datfile = data.frame(mat.or.vec(length(unique(Zfile1[,1]))-1,3))
print("head(datfile)1")
print(head(datfile))

print("dim(datfile)")
print(dim(datfile))

print("head(uniqid)")
print(head(uniqid))

colnames(datfile) = c("ID", "case", "peptide")
print("head(datfile)2")
print(head(datfile))



datfile$ID = uniqid
print("head(datfile)3")
print(head(datfile))

print("datfile$ID")
print(datfile$ID)

for(i in c(1:nrow(datfile))){
	datfile$case[i] = meta$group[which(meta$subject== datfile$ID[i])[1]]
}



print("head(Zfile1[,1])")
print(head(Zfile1[,1]))
print("head(datfile$ID)")
print(head(datfile$ID))

#	[1] "head(Zfile1[,1])"
#	[1] "id"         "4199"       "4199dup"    "1439301"    "1439301dup"
#	[6] "21129"     
#	[1] "head(datfile$ID)"
#	[1] "4199"    "1439301" "21129"   "1439701" "4465"    "1425001"






# Result File
pvalues = data.frame(mat.or.vec(ncol(Zfile1)-1, 5))
colnames(pvalues) = c("peptide", "species", "freq_case", "freq_control", "pval")


print("This loop takes quite a while. Haven't attempted to speed it up.")
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
		if(sum(dfl$peptide)==0 ) {
			pvalues$freq_case[pep_index] = 0
			pvalues$freq_control[pep_index] = 0
			pvalues$pval[pep_index] = NA
		} else {

			# If its present in everyone, report that
			if(sum(dfl$peptide)==nrow(dfl) ) {
				pvalues$freq_case[pep_index] = 1
				pvalues$freq_control[pep_index] = 1
				pvalues$pval[pep_index] = NA
			} else {

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

colnames(pvalues) = c("peptide", "species", paste0("freq_", groups_to_compare[1]),
	paste0("freq_", groups_to_compare[2]), "pval")


outfile=paste0(opt$output_dir, "/",
	gsub(" ","_",paste("Tile_Comparison",
		fs::path_ext_remove(basename(opt$zfilename)),
		"type",opt$type,
		paste(groups_to_compare[1:2],collapse="-"),
		"Z", Z,
		"sex",opt$sex, sep="-")), ".csv")


print(paste0("Writing ",outfile))
write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],
	outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


