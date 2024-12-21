#!/usr/bin/env Rscript


#	Creates manhattan plot pdf file for a given virus for all groups specified (utilizing the "type" column in the manifest file)
#	User provides name of virus, number of plots you want to have per page, ( I have it as 5 right now), and the groups you want to include (type column)
#	Indicates Public epitopes (as read from the "All.public_epitope_annotations.Zscores.csv" file) with red triangles.
#	Draws lines at Z = 3.5 and Z=10, and recodes any super large Z scores as Z = 15 for visualization sake.
#	Creates a PDF "Manhattan_plots_[VIRUS]" with indication of the groups compared.

library("optparse")

option_list = list(
	make_option(c("-v", "--virus"), type="character", default=NULL,
		help="virus name", metavar="character"),
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

if (is.null(opt$virus)){
	print_help(opt_parser)
	stop("virus name required.\n", call.=FALSE)
}



library(ggplot2)
library(gridExtra)



# Input parameters

# this can take any length of groups from the metadata file "type" column, no need to limit to 2.
groups_to_compare=c("PF Patient", "Endemic Control" , "Non Endemic Control")

#	groups_to_compare=unlist(strsplit(opt$groups, split = ","))

n_plots_per_page = 5





# Read in the Z-score file  (wihtout transpose and remove the transpose line)

Zfile = read.csv(paste(opt$working_dir, "Zscores.t.csv", sep = "/"), sep = ",", header=FALSE)

Zfile = data.frame(t(Zfile))

print("Read in the metadata file")
meta = read.csv( opt$manifest, sep= ",", header = TRUE)

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
public_eps = read.csv(paste(opt$working_dir, "All.public_epitope_annotations.Zscores.csv", sep = "/"), header = TRUE, sep = ",")
public_ep_id = public_eps$id
rm(public_eps)




print("Keep only this virus")
viral_ids = species_id$id[which(species_id$Species==opt$virus)]
Zfile2 = Zfile1[,c(1, which(Zfile1[1,] %in% viral_ids))]
rm(Zfile1)
# Order the viral IDs numerically
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

viral_plotfile = paste0(opt$working_dir,"/","Manhattan_plots-",gsub(" ","_",opt$virus), "-", gsub(" ","_",paste(groups_to_compare, collapse = "-")), ".pdf")
pdf(viral_plotfile, width = 7, height=(2*n_plots_per_page), onefile = TRUE)

for(stat in cc){
	statids = uniqid[which(uniqid %in% meta$subject[which(meta$group== stat)])]
	plots = list()
	plot_counter = 1
	for(indiv in statids){
		plot_df$Zscore = 0
		for(i in c(1:nrow(plot_df))){
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

