#!/usr/bin/env Rscript


#	Creates manhattan plot pdf file for a given virus for all groups specified (utilizing the "type" column in the manifest file)
#	User provides name of virus, number of plots you want to have per page, ( I have it as 5 right now), and the groups you want to include (type column)
#	Indicates Public epitopes (as read from the "All.public_epitope_annotations.Zscores.csv" file) with red triangles.
#	Draws lines at Z = 3.5 and Z=10, and recodes any super large Z scores as Z = 15 for visualization sake.
#	Creates a PDF "Manhattan_plots_[VIRUS]" with indication of the groups compared.


library("argparse")
args=commandArgs()
scriptname=sub("--file=", "", args[grepl("--file=", args)])
parser <- ArgumentParser(description=scriptname)
parser$add_argument("-g", "--groups_to_compare", type="character", required=TRUE, action="append",
	help="group to compare (use multiple times for each)", metavar="group")
parser$add_argument("-v", "--virus", type="character", default=NULL, required=TRUE,
	help="virus name", metavar="virus species")
parser$add_argument("-m", "--manifest", type="character", default=NULL, required=TRUE,
	help="manifest file name", metavar="manifest")
parser$add_argument("--output_dir", type="character", default="./",
	help="output dir [default=%(default)s]", metavar="directory")
parser$add_argument("--zfilename", type="character", default="out/Zscores.csv",
	help="zfilename [default=%(default)s]", metavar="Zscores file")
parser$add_argument("--public_eps_filename", type="character", default="out/All.public_epitope_annotations.Zscores.csv",
	help="public_eps_filename [default=%(default)s]", metavar="Public Eps file")
#parser$add_argument("-d", "--working_dir", type="character", default="./",
#	help="working dir [default=%(default)s]", metavar="directory")
opt <- parser$parse_args()


library(ggplot2)
library(gridExtra)
library(data.table)


# Input parameters
# this can take any length of groups from the metadata file "type" column, no need to limit to 2.

#groups_to_compare=opt$groups_to_compare
groups_to_compare=unlist(strsplit(opt$groups_to_compare, split = ","))
print("Comparing these groups")
print(groups_to_compare)

n_plots_per_page = 5


# Read in the Z-score file  (wihtout transpose and remove the transpose line)

#Zfile = read.csv(paste(opt$working_dir, "Zscores.t.csv", sep = "/"), sep = ",", header=FALSE)
#Zfile = data.frame(t(Zfile))
#Zfile = read.csv(opt$zfilename, sep = ",", header=FALSE)
Zfile <- data.frame(data.table::fread(opt$zfilename, sep = ",", header=FALSE))
print("head(Zfile)")
print(Zfile[1:5,1:5])

print("Read in the metadata file")
#meta = read.csv( opt$manifest, sep= ",", header = TRUE)
meta <- data.frame(data.table::fread(opt$manifest, sep = ",", header=TRUE))

# ## Code to create an id_species file so Jake doesn't need to append it each time
# mydat = Zfile[c(1,2),-c(1:3)]
# id_species = data.frame((t(mydat)))
# id_species = id_species[order(id_species$V1, decreasing = FALSE),]
# colnames(id_species) = c("id", "species")
# write.table(id_species, paste(mwd, "/ID_species.csv", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")


# If in the format of subject, type species, remove subject and type, and remove second row.
if("subject" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "subject")
	Zfile = Zfile[,-to_remove]
}

if("type" %in% Zfile[2,c(1:3)]){
	to_remove = which(Zfile[2,c(1:3)]== "type")
	Zfile = Zfile[,-to_remove]
}


print(Zfile[1:5,1:5])

print("Extract the peptide information")

#read.csv(paste(mwd, "/ID_species.csv", sep = ""), sep = ",", header = TRUE)
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


# This simply needs to be a list of peptides to work. Just need to know which peptides are deemed as "Public Epitopes"
#public_eps = read.csv(paste(opt$working_dir, "All.public_epitope_annotations.Zscores.csv", sep = "/"), header = TRUE, sep = ",")
#public_eps = read.csv(opt$public_eps_filename, header = TRUE, sep = ",")
public_eps <- data.frame(data.table::fread(opt$public_eps_filename, sep = ",", header=TRUE))
public_ep_id = public_eps$id
rm(public_eps)


print("Keep only this virus")
print("IT IS IMPORTANT THAT THE species column name is 'species' NOT 'Species'")
print(opt$virus)
viral_ids = species_id$id[which(species_id$species==opt$virus)]
print("viral_ids[1:4]")
print(viral_ids[1:4])
Zfile2 = Zfile1[,c(1, which(Zfile1[1,] %in% viral_ids))]
print("Zfile1[1:4,1:4]")
print(Zfile1[1:4,1:4])
rm(Zfile1)

print("Order the viral IDs numerically")
print(Zfile2[1:4,1:4])
#[1] "id"          "14078-01"    "14078-01dup" "14118-01"    "14118-01dup"
IDs = as.numeric(Zfile2[1, -c(1)])
IDs_sort = sort(IDs, decreasing = FALSE)
Zfile3 = Zfile2[,c(1,order(IDs,decreasing = FALSE)+1 )]
rm(Zfile2)

cc = groups_to_compare
plot_df = data.frame(mat.or.vec(length(IDs),4))
colnames(plot_df) = c("index", "ID", "Public", "Zscore")
plot_df$ID = IDs_sort
plot_df$index = c(1:nrow(plot_df))
for(i in c(1:nrow(plot_df))){
	ID = plot_df$ID[i]
	plot_df$Public[i] = ifelse(ID %in% public_ep_id, 1, 0)
}

#	opt$working_dir,"/",
print("Opening plot file")
viral_plotfile = paste0(
	opt$output_dir,"/",
	gsub(" ","_", paste("Manhattan_plots", opt$virus, paste(groups_to_compare, collapse = "-"),sep="-")),
	".pdf")

pdf(viral_plotfile, width = 7, height=(2*n_plots_per_page), onefile = TRUE)

for(stat in cc){
	print(paste0("for(stat in cc){ - ",stat))
	statids = uniqid[which(uniqid %in% meta$subject[which(meta$group== stat)])]
	plots = list()
	plot_counter = 1
	before1=Sys.time()
	for(indiv in statids){
		print(paste0("for(indiv in statids){ - ",indiv))
		plot_df$Zscore = 0
		print(paste0("Time : ",format(Sys.time(),"%Y%m%d%H%M%S")))
		before2=Sys.time()
		for(i in c(1:nrow(plot_df))){
			print(paste0("for(i in c(1:nrow(plot_df))){ - ",i))
			ID = plot_df$ID[i]
			Z = min(as.numeric(Zfile3[grep(indiv,Zfile3[,1]), which(Zfile3[1,]== ID)] ))
			if(is.na(Z)){
				Z = 0
			}
			if(is.infinite(Z)){
				Z = 15
			}
			if ( Z <0){
				Z = 0
			}
			if( Z> 15){
				Z = 15
			}

			# Implement something for the INFs

			plot_df$Zscore[i] = Z
		} # Loop over epitope IDs
		print("Done with the zscore min loop")
		print(paste0("Time : ",format(Sys.time(),"%Y%m%d%H%M%S")))
		after2=Sys.time()
		print(paste0(format(before2,"%Y%m%d%H%M%S")," - ",format(after2,"%Y%m%d%H%M%S")))
		print(difftime(after2,before2,units='min'))

		plot_df$Public = factor(plot_df$Public, levels = c("0", "1"))
		plots[[plot_counter]] = ggplot(plot_df, aes(x=index, y=Zscore, color=Public, shape = Public, size =Public)) +
			geom_point()+
			scale_color_manual(values=c('grey80','red')) +
			scale_size_manual(values = c(2, 4))+
			theme_classic()+
			theme(legend.position = "none")  +
			labs(title=paste0(indiv, " : ", stat),
				x="", y = "Minimum Z-score") +
			geom_hline(yintercept = (3.5), linetype="dashed", color = "blue") +ylim(c(0,17)) +
			geom_hline(yintercept = (10), linetype="dashed", color = "gold")
			plot_counter = plot_counter +1

	} # Loop over individuals in the group
	after1=Sys.time()
	print(difftime(after1,before1,units='min'))

	# This is an attempt at allowing an arbitrary number of plots be drawn, separating pages into only plots of max n_plots_per_page.
	# Not sure if the while logic may fail in some case I havent thought of yet.

	pcount = 1
	pmax = min(length(plots), (pcount + n_plots_per_page - 1))
	while(pcount <= pmax){
		myplot = grid.arrange(grobs = lapply(pcount:pmax, function(i) plots[[i]]), ncol=1, nrow=n_plots_per_page, top = opt$virus)
		print(myplot)
		pcount = pcount+n_plots_per_page
		pmax = min(length(plots), (pcount + n_plots_per_page - 1))
	}

} # Loop over Case Control

dev.off()

print(paste0("Writing ", viral_plotfile ))

