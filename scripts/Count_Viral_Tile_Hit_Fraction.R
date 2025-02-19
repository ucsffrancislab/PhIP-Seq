#!/usr/bin/env Rscript

#	merged with Case_Control_Seropositivity_Frac.R as it requires the output of this script.

#	Creates a table "Viral_Frac_Hits_Z_*.csv" indicating the fraction of tiles "hit" for each virus, for each sample, for the specified Z. (number tiles hit for viral species/total tiles associated with viral species)
#	This is a prerequisite to running the below Case_Control_Seropositivity_Frac.R script.
#	This table is computed for all samples on the plate, so no need to specify anything, if they are in the manifest, their viral fractions are computed.
#	We are using this ultimately as an alternate way to call a virus present, by thresholding the proportion of tiles present (see next file). This of course comes with many caveats, like tile uniqueness/homology.


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


print("Read in the metadata file")

#meta = read.csv(opt$manifest, sep= ",", header = TRUE)
meta <- data.frame(data.table::fread(opt$manifest, sep = ",", header=TRUE))

print("Read in the Z file")

#Zfile = read.csv(paste(opt$working_dir, "Zscores.t.csv", sep = "/"), sep = ",", header=FALSE)
#Zfile = read.csv(opt$zfilename, sep = ",", header=FALSE)
Zfile <- data.frame(data.table::fread(opt$zfilename, sep = ",", header=FALSE))
#Zfile = data.frame(t(Zfile))
print(Zfile[1:5,1:5])

#	[1] "Comparing these groups"
#	[1] "case"    "control"
#	[1] "Read in the metadata file"
#	[1] "Read in the Z file"

#	        V1           V2          V3                    V4                   V5
#	1        y            x          id                     1                   10
#	2  subject         type     species Papiine herpesvirus 2       Vaccinia virus
#	3 14078-01 glioma serum    14078-01   -0.1735488911468417  -0.1735488911468417
#	4 14078-01 glioma serum 14078-01dup   -0.2909274899485718  -0.2909274899485718
#	5 14118-01 glioma serum    14118-01  -0.20783649486361835 -0.20783649486361835

#	[1] "Extract the peptide information"
#	[1] "Unique samples to keep"
#	[1] "14078-01" "14118-01" "14127-01" "14142-01" "14206-01"
#	[1] "Unique species (column name is 'species' NOT 'Species')"
#	[1] "Papiine herpesvirus 2" "Vaccinia virus"        "Human herpesvirus 3"
#	[4] "Hepatitis B virus"     "Human herpesvirus 8"
#	[1] "Shell file for viral fractions"


# If in the format of subject, type species, remove subject and type, and remove second row.
if("subject" %in% Zfile[2,c(1:3)]){
	to_remove= which(Zfile[2,c(1:3)]== "subject")
	Zfile = Zfile[,-to_remove]
}
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

print("Zfile = Zfile[-2,]")
print(Zfile[1:5,1:5])

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

outfile=paste0(opt$output_dir, "/",
	gsub(" ","_", paste("Viral_Frac_Hits",
		fs::path_ext_remove(basename(opt$zfilename)),
		"Z", Z, paste(groups_to_compare[1:2],collapse="-"), sep="-")), ".csv")

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


outfile=paste0(opt$output_dir, "/",
	gsub(" ","_", paste("Viral_Sero_test_results",
		fs::path_ext_remove(basename(opt$zfilename)),
		paste(groups_to_compare[1:2],collapse="-"),
		"Vir_hit_frac", Vir_frac, "Z", Z, sep="-")), ".csv")

print(paste0("Writing ",outfile))
write.table(opvalues, outfile, col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)


