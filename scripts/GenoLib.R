
library(data.table)

# Read in the multiple manifest files.

read_multiple_manifests <- function(plates) {

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



