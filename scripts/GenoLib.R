
library(data.table)



select_subjects <- function(manifest,opt) {

	print("Selecting subjects")

	uniq_sub = unique(manifest$subject[which(manifest$group %in% c(opt$group1,opt$group2))])

	if ( opt$sex == "" ){
		print("Sex is not set so not filtering on sex.")
	} else {
		print(paste0("Sex is set to ",opt$sex,". Filtering"))
		uniq_sub = intersect(uniq_sub,unique(manifest$subject[which( manifest$sex==opt$sex )]))
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

