#!/usr/bin/env Rscript

#	Creates a table "Viral_Frac_Hits_Z_*.csv" indicating the fraction of tiles "hit" for each virus, for each sample, for the specified Z. (number tiles hit for viral species/total tiles associated with viral species)
#	This is a prerequisite to running the below Case_Control_Seropositivity_Frac.R script.
#	This table is computed for all samples on the plate, so no need to specify anything, if they are in the manifest, their viral fractions are computed.
#	We are using this ultimately as an alternate way to call a virus present, by thresholding the proportion of tiles present (see next file). This of course comes with many caveats, like tile uniqueness/homology.



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



# Compute fraction of virus specific tiles hit for each sample

# For each virus, calculate the number of tiles called positive (both reps Z > Z threshold), and divide by the total number of represented tiles for that virus.


# Input parameters



Z = 3.5

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


