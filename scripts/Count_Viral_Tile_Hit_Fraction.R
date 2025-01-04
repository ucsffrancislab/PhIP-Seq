#!/usr/bin/env Rscript

#	merged with Case_Control_Seropositivity_Frac.R as it requires the output of this script.

#	Creates a table "Viral_Frac_Hits_Z_*.csv" indicating the fraction of tiles "hit" for each virus, for each sample, for the specified Z. (number tiles hit for viral species/total tiles associated with viral species)
#	This is a prerequisite to running the below Case_Control_Seropositivity_Frac.R script.
#	This table is computed for all samples on the plate, so no need to specify anything, if they are in the manifest, their viral fractions are computed.
#	We are using this ultimately as an alternate way to call a virus present, by thresholding the proportion of tiles present (see next file). This of course comes with many caveats, like tile uniqueness/homology.


library("optparse")

option_list = list(
  make_option(c("-z", "--zscore"), type="double", default=3.5,
    help="Zscore threshold", metavar="character"),
	make_option(c("-a", "--group1"), type="character", default=NULL,
		help="First group to compare", metavar="character"),
	make_option(c("-b", "--group2"), type="character", default=NULL,
		help="Second group to compare", metavar="character"),
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

groups_to_compare = c(opt$group1,opt$group2)

print("Comparing these groups")
print(groups_to_compare)







# Compute fraction of virus specific tiles hit for each sample

# For each virus, calculate the number of tiles called positive (both reps Z > Z threshold), and divide by the total number of represented tiles for that virus.


# Input parameters



#Z = 3.5
Z = opt$zscore

print("Read in the metadata file")

meta = read.csv(opt$manifest, sep= ",", header = TRUE)

print("Read in the Z file")

Zfile = read.csv(paste(opt$working_dir, "Zscores.t.csv", sep = "/"), sep = ",", header=FALSE)
Zfile = data.frame(t(Zfile))


# If in the format of subject, type species, remove subject and type, and remove second row.
if("subject" %in% Zfile[2,c(1:3)]){
	to_remove= which(Zfile[2,c(1:3)]== "subject")
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
uniqid = unique(meta$subject)
print(uniqid[1:5])

to_keep = 1
for(u in uniqid){
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
virfracs = data.frame(mat.or.vec(length(uniqid), length(unique(species_id$species))+1))
colnames(virfracs) = c("id", unique(species_id$species))
virfracs$id = uniqid

for( j in c(2:ncol(virfracs))){
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

outfile=paste0(opt$working_dir, "/", "Viral_Frac_Hits_Z_", Z, ".csv")
print(paste0("Writing ",outfile))
write.table(virfracs, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

#	warnings()
#	1: In FUN(newX[, i], ...) : no non-missing arguments to min; returning Inf
#	2: In FUN(newX[, i], ...) : no non-missing arguments to min; returning Inf




#	merged 2 scripts since one needs the output of the other



#	Assesses differences in the proportion of samples with virus called present between two user specified groups (from the manifest "type" column).
#	Uses our own definition of seropositivity (Virus present if >"Vir_frac" proportion of tiles are hit at Z > "Z_thresh"). I typically have used Vir_frac = 0.02 or 0.05 based on eyeballing the Viral_Frac_Hits table generated above.
#	Outputs a file in the same test directory called "Viral_Sero_test_results*" with indicators of the groups and parameters used.



Vir_frac = 0.05
#virfracfilename = paste0("Viral_Frac_Hits_Z_",Z,".csv")

print("Read in the metadata file")

meta = read.csv( opt$manifest, sep= ",", header = TRUE)

#	don't need to read it as still in memory from above.
#	DIFFERENT VARIABLE NAME
#print("Read in the VirFrac file")
#vir_fracs = read.csv(paste(opt$working_dir, virfracfilename, sep = "/"), header = TRUE)

# # Read in the viral score file
# vs = read.csv(paste(mwd, "/", virfilename, sep =""),sep = ",", header= FALSE )
# vir_score = data.frame(t(vs))
# colnames(vir_score) = vir_score[1,]
# vir_score = vir_score[-1,]


print("Unique samples to keep")
uniqid = unique(meta$subject[which(meta$group %in% groups_to_compare)])
#vir_score = vir_score[which(vir_score$id %in% uniqid),]
virfracs = virfracs[which(virfracs$id %in% uniqid),]

cases = unique(meta$subject[which(meta$group %in% groups_to_compare[1])])
controls = unique(meta$subject[which(meta$group %in% groups_to_compare[2])])

print("Create a shell file for analysis")
pvalues = data.frame(mat.or.vec(ncol(virfracs)-1, 4))
colnames(pvalues) = c( "species", "freq_case", "freq_control", "pval")
pvalues$species = colnames(virfracs)[-1]


for(i in c(1:nrow(pvalues))){
	species = pvalues$species[i]
	n_cases = length(cases)
	n_control = length(controls)

	# For each case, determine number of ids that are >3.5 in both reps

	n_case_success = length(which(as.numeric(virfracs[which(virfracs$id %in% cases), which(colnames(virfracs) == species)]) > Vir_frac))
	n_control_success = length(which(as.numeric(virfracs[which(virfracs$id %in% controls), which(colnames(virfracs) == species)]) > Vir_frac))

	prop = prop.test(c(n_case_success, n_control_success), c(n_cases, n_control), p = NULL, alternative = "two.sided", correct = TRUE)
	pvalues$freq_case[i] = n_case_success/n_cases
	pvalues$freq_control[i] = n_control_success/n_control
	pvalues$pval[i] = prop$p.value
}

colnames(pvalues) = c( "species", paste0("freq_", groups_to_compare[1]), paste0("freq_", groups_to_compare[2]), "pval")
opvalues = pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),]


outfile=paste0(opt$working_dir, "/",
	gsub(" ","_", paste("Viral_Sero_test_results", paste(groups_to_compare[1:2],collapse="-"),
		"Vir_hit_frac", Vir_frac, "Z", Z, sep="-")), ".csv")

print(paste0("Writing ",outfile))
write.table(opvalues, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


