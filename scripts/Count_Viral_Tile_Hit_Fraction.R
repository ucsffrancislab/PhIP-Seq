#!/usr/bin/env Rscript

#	merged with Case_Control_Seropositivity_Frac.R as it requires the output of this script.

#	Creates a table "Viral_Frac_Hits_Z_*.csv" indicating the fraction of tiles "hit" for each virus,
#	for each sample, for the specified Z.
#	(number tiles hit for viral species/total tiles associated with viral species)
#	This is a prerequisite to running the below Case_Control_Seropositivity_Frac.R script.
#	This table is computed for all samples on the plate, so no need to specify anything,
#		if they are in the manifest, their viral fractions are computed.
#	We are using this ultimately as an alternate way to call a virus present,
#		by thresholding the proportion of tiles present (see next file).
#	This of course comes with many caveats, like tile uniqueness/homology.

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
parser$add_argument("--zfilename", type="character", default="out/Zscores.csv",
	help="zfilename [default=%(default)s]", metavar="Zscores file")
opt <- parser$parse_args()

dir.create(opt$output_dir,showWarnings=F)

# For each virus, make a call of positive or negative based on at least 5% of possible tiles hitting, then measure proportion.
# Input parameters

groups_to_compare = c(opt$group1,opt$group2)

print("Comparing these groups")
print(groups_to_compare)


# Compute fraction of virus specific tiles hit for each sample

# For each virus, calculate the number of tiles called positive (both reps Z > Z threshold), and divide by the total number of represented tiles for that virus.

# Input parameters

Z = opt$zscore

library(data.table)


print("Read in the manifest file")
manifest <- data.frame(data.table::fread(opt$manifest, sep = ",", header=TRUE))


results=read_zfile(opt$zfilename)
species_id = results$species
Zfile = results$zfile
rm(results)


print("Unique samples to keep")
uniq_sub = select_subjects(manifest,opt)


to_keep = 1
for(u in uniq_sub){
	possible_ids = grep(u, Zfile[,1])
	mids = Zfile[possible_ids,1]
	locs = grepl("dup", mids)
	to_keep = c(to_keep,possible_ids[c(which(locs==FALSE)[1],which(locs==TRUE)[1]) ] )
}
Zfile1 = Zfile[to_keep,]
rm(Zfile)

print("Unique species (column name is 'species' NOT 'Species')")
uniq_spec = unique(species_id$species)
print(uniq_spec[1:5])

print("Shell file for viral fractions")
virfracs = data.frame(mat.or.vec(length(uniq_sub), length(unique(species_id$species))+1))
colnames(virfracs) = c("id", unique(species_id$species))
virfracs$id = uniq_sub

print(paste("Looping",ncol(virfracs)))
for( j in c(2:ncol(virfracs))){
	print(paste("Looping",j,":",ncol(virfracs)))
	sp = colnames(virfracs)[j]
	vir_ids = species_id$id[which(species_id$species == sp)]
	sfile = Zfile1[, c(1, which(Zfile1[1,] %in% vir_ids))]

	for(i in c(1:nrow(virfracs))){
		vf = sfile[grep(virfracs$id[i], sfile[,1]), ]
		if(ncol(vf)>2){
			n_positive = length(which( as.numeric(apply((vf[,-1]),2,min)) >Z))
		}else{
			n_positive = ifelse(as.numeric(min(vf[,2])) >Z, 1, 0)
		}
		virfracs[i, j] = n_positive/length(vir_ids)
	}
}

outfile=paste0(opt$output_dir, "/",
	gsub(" ","_", paste("Viral_Frac_Hits",
		fs::path_ext_remove(basename(opt$zfilename)),
		"type",opt$type,
		paste(groups_to_compare[1:2],collapse="-"),
		"Z", Z,
		"sex",opt$sex, sep="-")), ".csv")

print(paste0("Writing ",outfile))
write.table(virfracs, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

#	warnings()
#	1: In FUN(newX[, i], ...) : no non-missing arguments to min; returning Inf
#	2: In FUN(newX[, i], ...) : no non-missing arguments to min; returning Inf


#	merged 2 scripts since one needs the output of the other


#	Assesses differences in the proportion of samples with virus called present between two user
#	specified groups (from the manifest "type" column).
#	Uses our own definition of seropositivity (Virus present if >"Vir_frac" proportion of tiles are hit at Z > "Z_thresh").
#	I typically have used Vir_frac = 0.02 or 0.05 based on eyeballing the Viral_Frac_Hits table generated above.
#	Outputs a file in the same test directory called "Viral_Sero*" with indicators of the groups and parameters used.


Vir_frac = 0.05


virfracs = virfracs[which(virfracs$id %in% uniq_sub),]
print(paste("length(virfracs) :",length(virfracs)))

cases = unique(manifest$subject[which(manifest$group %in% groups_to_compare[1])])
cases = intersect(cases,uniq_sub)
print(paste("length(cases) :",length(cases)))
controls = unique(manifest$subject[which(manifest$group %in% groups_to_compare[2])])
controls = intersect(controls,uniq_sub)
print(paste("length(controls) :",length(controls)))

print("Create a shell file for analysis")
pvalues = data.frame(mat.or.vec(ncol(virfracs)-1, 4))
colnames(pvalues) = c("species", "freq_case", "freq_control", "pval")
pvalues$species = colnames(virfracs)[-1]


for(i in c(1:nrow(pvalues))){
	print(paste("Looping:",i,":",nrow(pvalues)))
	species = pvalues$species[i]
	n_cases = length(cases)
	n_control = length(controls)

	# For each case, determine number of ids that are >3.5 in both reps

	n_case_success = length(which(as.numeric(virfracs[which(virfracs$id %in% cases),
		which(colnames(virfracs) == species)]) > Vir_frac))
	n_control_success = length(which(as.numeric(virfracs[which(virfracs$id %in% controls),
		which(colnames(virfracs) == species)]) > Vir_frac))

	prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control), p = NULL, alternative = "two.sided", correct = TRUE)
	pvalues$freq_case[i] = n_case_success/n_cases
	pvalues$freq_control[i] = n_control_success/n_control
	pvalues$pval[i] = prop$p.value
}

colnames(pvalues) = c( "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]


outfile=paste0(opt$output_dir, "/",
	gsub(" ","_", paste("Viral_Sero",
		fs::path_ext_remove(basename(opt$zfilename)),
		"type",opt$type,
		paste(groups_to_compare[1:2],collapse="-"),
		"Z", Z, 
		"sex",opt$sex,
		"Vir_hit_frac", Vir_frac, sep="-")), ".csv")

print(paste0("Writing ",outfile))
write.table(opvalues, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


