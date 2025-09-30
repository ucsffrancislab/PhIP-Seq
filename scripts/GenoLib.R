
library(data.table)


build_datfile = function(uniq_sub,opt){

	print("Building datfile")

	#datfile = data.frame(mat.or.vec(length(uniq_sub),6))
	#colnames(datfile) = c("ID", "case", "peptide", "sex", "age", "plate")
	datfile = data.frame(mat.or.vec(length(uniq_sub),5))
	colnames(datfile) = c("ID", "case", "sex", "age", "plate")
	datfile$ID = uniq_sub
	for(i in c(1:nrow(datfile))){
		print(paste("Looping:",i,":",nrow(datfile)))
		man_loc = which(manifest$subject== datfile$ID[i])[1]
		#datfile$case[i] = ifelse(manifest$group[man_loc] == groups_to_compare[1], 1, 0)
		datfile$case[i] = ifelse(manifest$group[man_loc] == opt$group1, 1, 0)
		datfile$age[i] = manifest$age[man_loc]
		datfile$sex[i] = manifest$sex[man_loc]
		datfile$plate[i] = manifest$plate[man_loc]
	}
	datfile$age = as.numeric(datfile$age)
	datfile$sex = as.factor(datfile$sex)
	datfile$plate = as.factor(datfile$plate)

	return(datfile)
}

read_zfile = function(zfilename) {

	print(paste("Read in the Z file",zfilename))

	Zfile <- data.frame(data.table::fread(zfilename, sep = ",", header=FALSE))
	print(Zfile[1:5,1:4])
	#       V1                 V2        V3                  V4
	#1       y                  x        id                  10
	#2 subject               type   species      Vaccinia virus
	#3  108741 ALL maternal serum    108741 -0.3633046359217093
	#4  108741 ALL maternal serum 108741dup -0.3669457592878151
	#5   14922       glioma serum     14922 -0.2185790682273388
	
	# If in the format of subject, type species, remove subject and type, and remove second row.
	if("subject" %in% Zfile[2,c(1:3)]){
		to_remove= which(Zfile[2,c(1:3)]== "subject")
		Zfile = Zfile[,-to_remove]
	}
	if("type" %in% Zfile[2,c(1:3)]){
		to_remove = which(Zfile[2,c(1:3)]== "type")
		Zfile = Zfile[,-to_remove]
	}
	
	print(Zfile[1:5,1:4])
	#         V3                  V4                   V5                  V6
	#1        id                  10                  100                1000
	#2   species      Vaccinia virus  Human herpesvirus 3   Hepatitis B virus
	#3    108741 -0.3633046359217093  -0.5834907519966614 -0.3633046359217093
	#4 108741dup -0.3669457592878151  -0.5919839162855073 -0.3669457592878151
	#5     14922 -0.2185790682273388 -0.40674341138804637 -0.2185790682273388	
	
	species_id = data.frame(t(Zfile[c(1:2),]))
	print(species_id[1:3,1:2])
	#      X1                  X2
	#V3    id             species
	#V4    10      Vaccinia virus
	#V5   100 Human herpesvirus 3

	colnames(species_id) = species_id[1,]
	print(species_id[1:3,1:2])
	#      id             species
	#V3    id             species
	#V4    10      Vaccinia virus
	#V5   100 Human herpesvirus 3

	species_id = species_id[-1,]
	print(species_id[1:3,1:2])
	#      id             species
	#V4    10      Vaccinia virus
	#V5   100 Human herpesvirus 3
	#V6  1000   Hepatitis B virus

	Zfile = Zfile[-2,]
	print(Zfile[1:5,1:4])
	#         V3                  V4                   V5                  V6
	#1        id                  10                  100                1000
	#3    108741 -0.3633046359217093  -0.5834907519966614 -0.3633046359217093
	#4 108741dup -0.3669457592878151  -0.5919839162855073 -0.3669457592878151
	#5     14922 -0.2185790682273388 -0.40674341138804637 -0.2185790682273388
	#6  14922dup -0.2617924942197322  -0.4092362704232559 -0.2617924942197322

	colnames(Zfile) = Zfile[1,]
	print(Zfile[1:5,1:4])

	return(list(zfile = Zfile, species = species_id))

}


log_reg = function(df,logitmodel,label){

	# A simple model that simply adjusts for plate/batch in the model. When the number
	#	of plates becomes large, a mixed effects regression model should be considered.
	# as there are likely differences in the peptide calling sensitivity between plates,
	#	and so peptide probably has different associations with case based on plate.

	logit_fun = glm(as.formula(logitmodel), data = df, family=binomial(link="logit"))
	go = summary(logit_fun)

	cat(capture.output(print(go)), file = logname, append = TRUE, sep = "\n")

	#beta = go$coefficients[2,1]
	#se = go$coefficients[2,2]
	#pval = go$coefficients[2,4]
	#	sometimes "Coefficients: (1 not defined because of singularities)"
	#	and peptide doesn't exist. This would then return the age coefficients.
	beta <- if(label %in% rownames(go$coefficients)) go$coefficients[label,'Estimate']   else NA
	se <-   if(label %in% rownames(go$coefficients)) go$coefficients[label,'Std. Error'] else NA
	pval <- if(label %in% rownames(go$coefficients)) go$coefficients[label,'Pr(>|z|)']   else NA

	return(c(beta, se, pval))
}



select_subjects <- function(manifest,opt) {

	print("Selecting subjects")

	uniq_sub = unique(manifest$subject[which(manifest$group %in% c(opt$group1,opt$group2))])

	if ( opt$sex == "" ){
		print("Sex is not set so not filtering on sex.")
	} else {
		print(paste0("Sex is set to ",opt$sex,". Filtering"))
		uniq_sub = intersect(uniq_sub,unique(manifest$subject[which( manifest$sex==opt$sex )]))
	}

#	if ( opt$study == "" ){
#		print("Study is not set so not filtering on study.")
#	} else {
#		print(paste0("Study is set to ",opt$study,". Filtering"))
#		uniq_sub = intersect(uniq_sub,unique(manifest$subject[which( manifest$study==opt$study )]))
#	}

	if ( length(opt$study) == 0 ){
		print("Study is not set so not filtering on study.")
	} else {
		print(paste0("Study is set to ",paste(opt$study),". Filtering"))
		uniq_sub = intersect(uniq_sub,unique(manifest$subject[which( manifest$study %in% opt$study )]))
	}

	if ( opt$type == "" ){
		print("Type is not set so not filtering on type.")
	} else {
		print(paste0("Type is set to ",opt$type,". Filtering"))
		uniq_sub = intersect(uniq_sub,unique(manifest$subject[which( manifest$type==opt$type )]))
	}

	print(paste("Length select subjects:", length(uniq_sub) ) )
	print(paste("Unique subjects:", paste(uniq_sub, collapse=", ") ))

	return(uniq_sub)
}



# Read in the multiple manifest files.

read_multiple_manifests <- function(plates) {
	print("Reading multiple manifests")
	print(plates)

	mfs = list()

	for(i in c(1:length(plates))){

		# Find the manifest file for the given plate. Requires only ONE manifest file per plate folder
		mfname = list.files(plates[i], pattern="manifest", full.names=TRUE)
		if(length(mfname)!=1){
			print(paste0(plates[i], " needs a single manifest file!"))
		}

		# read in the manifest file
		mf <- data.frame(data.table::fread(mfname, sep = ",", header=TRUE)) # 50x faster

		# Create a categorical variable, assign all of these the same number to indicate plate.

		# The plate number is now a column in the manifest file
		# Later in this script it is used in a formula and can't all the be the same
		# However, it is used as the index and not the actual plate number.
		mf$plate = i 
		mfs[[i]] = mf
	}

	# Create an aggregate metadata file. This requires identical column structure in the files.

	manifest = Reduce(rbind, mfs)
	# can get rid of mfs list.
	rm(mfs)

	return(manifest)
}



print("Loaded GenoLib")

