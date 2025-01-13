#!/usr/bin/env Rscript


#	The multiplate script to do analysis for our own viral calls using a specified fraction hit.
#	Multi_Plate_Case_Control_VirHitFrac_Seropositivity_Regression.R
#	Uses a logistic regression framework to assess the association between the presence of each virus and user specified case control status.
#	Uses our call of viral positivity using the viral hit frac file. (Which you need to run in each plate prior to running this) Z_thresh and Vir_frac are necessary inputs, Z to know to pull the right file, and Vir_frac to make 1/0 viral calls.
#	Takes as input a list of directories (e.g. plates), to use in the analysis, the directories should contain the Viral_Frac_Hits_Z-*.csv and the manifest* files.
#	Only viral names found in all plates (just in the column names, not requirement of actual presence) are included in the analysis
#	The logistic regression adjusts for age (continuous), sex (M/F factor), and plate (1..n plates, factor)
#	The output (to a user specified directory owd , is a list of viral species, ordered by ascending p-value. It includes species, frequency of virus in the first group, frequency of the virus in the second group, and the beta, standard error, and p-value from the multivariable logistic regression (NA if there is no variation in the presence of the virus overall, e.g. all 0 or all 1).
#	The file will also output a logfile with similar naming convention, which includes the details of the plates used in the analysis, sample sizes, etc.



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
parser$add_argument("-p", "--plate", type="character", required=TRUE, action="append",
	help="plate to compare (use multiple times for each)", metavar="group")
parser$add_argument("-o", "--output_dir", type="character", default="./",
	help="output dir [default=%(default)s]", metavar="directory")
opt <- parser$parse_args()


groups_to_compare = c(opt$group1,opt$group2)
print("Comparing these groups")
print(groups_to_compare)

plates=opt$plate
print("Comparing these plates")
print(plates)

owd=opt$output_dir
print("Output dir")
print(owd)
dir.create(owd,showWarnings=F)

# Multi plate Logistic regression for seropositivity (hit Frac) on case/control status, adjusting for age, sex, and plate


Z_thresh = opt$zscore
Vir_frac = 0.05

#	DO NOT INCLUDE THE DIRNAME
virfracfilename = paste0(
	gsub(" ","_", paste("Viral_Frac_Hits_Z", Z_thresh,
	paste(groups_to_compare[1:2],collapse="-"), sep="-")), ".csv")

date=format(Sys.Date(),"%Y%m%d")

output_base = paste0(owd, "/", gsub(" ","_",
	paste(date, "Multiplate_VirHitFrac_Seropositivity_Comparison",
	paste(groups_to_compare, collapse="-"),"test_results",Z_thresh, sep="-")))

# Log the parameter choices into a logfile
logname = paste0(output_base,'.log')

cat("Multi plate Logistic regression for presence of virus on case/control status, adjusting for age, sex, and plate. Using our calculated hit fraction for virus calling.",file=logname,sep="\n")

cat("\nPlates used in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(plates))){
  cat(plates[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)

cat("\nGroups compared in this analysis:", file = logname, append = TRUE, sep = "\n")
for(i in c(1:length(groups_to_compare))){
  cat(groups_to_compare[i], file = logname, append= TRUE, sep = "\n")
}
cat("\n", file = logname, append = TRUE)

virfiles = list()
mfs = list()
# Read in multiple plate seropositivity files.
for(i in c(1:length(plates))){
  virfile = read.csv(paste0(plates[i], "/", virfracfilename), header = TRUE, sep = ",")


  virfiles[[i]] = virfile


  mfname = list.files(plates[i], pattern="manifest", full.names=TRUE)
  if(length(mfname)!=1){
    print(paste0(plates[i], " needs a single manifest file!"))
  }

  # read in the manifest file
  mf = read.csv(paste(mfname, sep = ""), sep= ",", header = TRUE)
  # Create a categorical variable, assign all of these the same number to indicate plate.
  mf$plate = i
  mfs[[i]] = mf

}# close loop over plates.


manifest = Reduce(rbind, mfs)
# can get rid of mfs list.
rm(mfs)



# Identify the unique subjects to include in the analyses.

uniq_sub = unique(manifest$subject[which(manifest$group %in% groups_to_compare)])

cat(paste0("\nTotal number of included subjects: ", length(uniq_sub)),
	file = logname, append = TRUE, sep = "\n")

# Identify the viruses that are in both files, and subset each posfile to that, with maintained order.
common_virs = Reduce(intersect, lapply(virfiles,function(x) colnames(x)))

for(i in c(1:length(plates))){
	virfiles[[i]] = virfiles[[i]][,common_virs]
}
common_virs = common_virs[-c(1)]

cat(paste0("\nTotal number of included viruses: ", length(common_virs)),
	file = logname, append = TRUE, sep = "\n")


# Convert every virus call into a binary 0/1 call based on >1 hits (Jake already filtered for public epitopes).

cat("\nStart converting viral calls to  binary calls", file = logname, append = TRUE, sep = "\n")

viral_calls = data.frame(mat.or.vec(length(uniq_sub), length(common_virs)+1))
colnames(viral_calls) = c("ID", common_virs)
viral_calls$ID = uniq_sub

# Loop over every person, and convert all of their peptide Z scores into 0's and 1's.
for(i in c(1:nrow(viral_calls))){
	id = viral_calls$ID[i]
	# Get plate to pull from
	mp = manifest$plate[which(manifest$subject==id)][[1]]

	# Extract the rows containing the sample name
	rloc = which(virfiles[[mp]][,1] == id)[1]

	myCounts = data.frame(t(virfiles[[mp]][rloc, which(colnames(virfiles[[mp]]) %in% common_virs)]))
	myCounts[,1]=as.numeric(myCounts[,1])
	if(length(which((common_virs == row.names(myCounts))==FALSE))>0){
		cat(paste("\nError:", id, ":", "plate", mp, ":Viruses out of order, need to ensure they are ordered the same as the common_virs vector. This should never appear.", sep = " "), file=logname, append = TRUE, sep= "\n")
	}

	viral_calls[i, -1] = ifelse(myCounts > Vir_frac, 1, 0)
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(virfiles)


# Create a shell file for analysis, Leaves the virus column blank, we will repopulate this with every peptide.

datfile = data.frame(mat.or.vec(length(uniq_sub),6))
colnames(datfile) = c("ID", "case", "virus", "sex", "age", "plate")
datfile$ID = uniq_sub
datfile$virus = NA
for(i in c(1:nrow(datfile))){
	man_loc = which(manifest$subject== datfile$ID[i])[1]
	datfile$case[i] = ifelse(manifest$group[man_loc] == groups_to_compare[1], 1, 0)
	datfile$age[i] = manifest$age[man_loc]
	datfile$sex[i] = manifest$sex[man_loc]
	datfile$plate[i] = manifest$plate[man_loc]
}


datfile$age = as.numeric(datfile$age)
datfile$sex = as.factor(datfile$sex)
datfile$plate = as.factor(datfile$plate)



#----- Shell function for logistic regression analysis.
log_reg = function(df){

	# A simple model that simply adjusts for plate/batch in the model. When the number of plates becomes large, a mixed effects regression model should be considered.
	# as there are likely differences in the virus calling sensitivity between plates, and so virus probably has different associations with case based on plate.
	logitmodel = "case~ virus +age + sex + plate"

	logit_fun = glm(as.formula(logitmodel), data = df, family=binomial(link="logit"))

	go= summary(logit_fun)
	beta = go$coefficients[2,1]
	se = go$coefficients[2,2]
	pval = go$coefficients[2,4]
	return(c(beta, se, pval))
}
#-------


# Result File
pvalues = data.frame(mat.or.vec(length(common_virs), 6))
colnames(pvalues) = c( "species",  "freq_case", "freq_control", "beta", "se", "pval")
pvalues$species = common_virs

n_case = length(which(datfile$case==1))
n_control = length(which(datfile$case==0))

cat(paste("\nTotal number of ", groups_to_compare[1], ": ", n_case, sep =""), file = logname, append = TRUE, sep = "\n")
cat(paste("\nTotal number of ", groups_to_compare[2], ": ", n_control, sep =""), file = logname, append = TRUE, sep = "\n")


# Loop over viruses, populate the datfile, and run analyses.
cat("\nStart loop over virus logistic regression analysis:", file = logname, append = TRUE, sep = "\n")

for(i in c(1:length(common_virs))){

	# Extract all virus information into the new datfile.
	datfile$virus = viral_calls[, which(colnames(viral_calls)== common_virs[i])]

	# Calculate the frequency of the peptide in each group
	n_case_pos = length(which(datfile$virus==1 & datfile$case==1))
	n_control_pos = length(which(datfile$virus==1 & datfile$case==0))

	pvalues$freq_case[i] = n_case_pos/n_case
	pvalues$freq_control[i] = n_control_pos/n_control

	if( (n_case_pos +n_control_pos) %in% c(0, n_case+n_control)){
		pvalues$beta[i]= NA
		pvalues$se[i] = NA
		pvalues$pval[i] = NA
	}else{
		results =log_reg(datfile)
		pvalues$beta[i]= results[1]
		pvalues$se[i] = results[2]
		pvalues$pval[i] = results[3]
	}
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

colnames(pvalues) = c( "species",
	paste0("freq_", groups_to_compare[1]),
	paste0("freq_", groups_to_compare[2]),
	"beta", "se", "pval")

write.table(pvalues[order(pvalues$pval,decreasing = FALSE, na.last = TRUE),],
	paste0(output_base,'.csv'),
	col.names = TRUE, sep = ",", row.names=FALSE, quote= FALSE)

cat("\n Analysis complete.", file = logname, append = TRUE, sep = "\n")

