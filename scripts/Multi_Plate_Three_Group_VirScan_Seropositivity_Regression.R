#!/usr/bin/env Rscript

#   Three-group multinomial logistic regression for PhIP-seq viral seropositivity
#   using VirScan's calling method (public epitope + >1 unique read hit).
#
#   Extends Multi_Plate_Case_Control_VirScan_Seropositivity_Regression.R from
#   two groups to three by replacing binary logistic regression (glm / binomial)
#   with multinomial logistic regression (nnet::multinom).  The original two-group
#   scripts are unchanged; this script lives alongside them.
#
#   Model
#   -----
#   For every virus, fit:
#       group ~ virus [+ age] [+ sex] [+ plate]
#   where group is a three-level factor with group1 as the reference.
#   This yields two sets of coefficients per virus:
#       group2 vs. group1
#       group3 vs. group1
#
#   Output columns (per virus, sorted by ascending min p-value):
#       species,
#       freq_<group1>, freq_<group2>, freq_<group3>,
#       beta_<group2>_vs_<group1>, se_<group2>_vs_<group1>, pval_<group2>_vs_<group1>,
#       beta_<group3>_vs_<group1>, se_<group3>_vs_<group1>, pval_<group3>_vs_<group1>
#
#   Usage example:
#       module load r
#       Multi_Plate_Three_Group_VirScan_Seropositivity_Regression.R \
#           -z 5 --type "pemphigus serum" \
#           -a "PF Patient" -b "Non Endemic Control" -c "Healthy Control" \
#           --sfile_basename seropositive.5.csv \
#           -o /path/to/output \
#           -p /path/to/out.plate13 -p /path/to/out.plate14 \
#           -p /path/to/out.plate21
#
#   Dependencies: argparse, data.table, nnet, fs
#   Sourced libs:  GenoLib.R       (read_multiple_manifests)
#                  ThreeGroupLib.R (build_datfile_three_group, select_subjects_three_group,
#                                   multinom_reg, build_formula_three_group,
#                                   init_log_three_group, write_results_three_group)

options(warn = 1)

library("argparse")
args       <- commandArgs()
scriptname <- sub("--file=", "", args[grepl("--file=", args)])
scriptdir  <- dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(paste(scriptdir, "GenoLib.R",       sep = "/"))
source(paste(scriptdir, "ThreeGroupLib.R", sep = "/"))

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = scriptname)
parser$add_argument("--study", type = "character", action = "append",
    help = "limit study", metavar = "study")
parser$add_argument("-t", "--type", type = "character", default = "",
    help = "limit type", metavar = "type")
parser$add_argument("-a", "--group1", type = "character", required = TRUE,
    help = "first group (reference level)", metavar = "group")
parser$add_argument("-b", "--group2", type = "character", required = TRUE,
    help = "second group (contrasted against group1)", metavar = "group")
parser$add_argument("-c", "--group3", type = "character", required = TRUE,
    help = "third group (contrasted against group1)", metavar = "group")
parser$add_argument("-z", "--zscore", type = "double", default = 3.5,
    help = "Z-score threshold [default %(default)s]", metavar = "double")
parser$add_argument("-s", "--sex", type = "character", default = "",
    help = "limit sex", metavar = "sex")
parser$add_argument("-p", "--plate", type = "character", required = TRUE,
    action = "append",
    help = "plate directory (repeat for each plate)", metavar = "directory")
parser$add_argument("-o", "--output_dir", type = "character", default = "./",
    help = "output directory [default %(default)s]", metavar = "directory")
parser$add_argument("--sfile_basename", type = "character",
    default = "seropositive.csv",
    help = "seropositive file basename [default %(default)s]",
    metavar = "seropositive file basename")
parser$add_argument("--keep_all_ids", action = "store_true",
    help = "keep all sample IDs, not just those ending in '_B'")
parser$add_argument("--ignore_plate", action = "store_true",
    help = "exclude plate from the regression formula")
opt <- parser$parse_args()

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

groups_to_compare <- c(opt$group1, opt$group2, opt$group3)
print("Comparing these groups"); print(groups_to_compare)

plates <- opt$plate
print("Comparing these plates"); print(plates)

owd <- opt$output_dir
print("Output dir"); print(owd)
dir.create(owd, showWarnings = FALSE, recursive = TRUE)

library(data.table)

Z <- opt$zscore

output_base <- paste0(owd, "/", gsub(" ", "_",
    paste("Multiplate_Three_Group_VirScan_Seropositivity_Comparison",
          fs::path_ext_remove(basename(opt$sfile_basename)),
          paste(opt$study, collapse = "-"), opt$type,
          paste(groups_to_compare, collapse = "-"),
          "Z", Z,
          "sex", opt$sex, sep = "-")))

logname <- paste0(output_base, ".log")

init_log_three_group(
    logname  = logname,
    analysis = paste0(
        "Multi-plate three-group multinomial logistic regression for viral ",
        "seropositivity (VirScan method), adjusting for age, sex, and plate.\n",
        "Reference group (group1): ", opt$group1, "\n",
        "Z-score threshold = ", Z),
    plates   = plates,
    groups   = groups_to_compare
)

# ---------------------------------------------------------------------------
# Read per-plate VirScan seropositive files
# ---------------------------------------------------------------------------

posfiles <- list()

for (i in seq_along(plates)) {
    posfile <- data.frame(
        data.table::fread(paste(plates[i], opt$sfile_basename, sep = "/"),
                          sep = ",", header = TRUE),
        check.names = FALSE)

    posfiles[[i]] <- if (opt$keep_all_ids) posfile else posfile[grep("_B$", posfile$id), ]
    rm(posfile)
}

manifest <- read_multiple_manifests(plates)
uniq_sub <- select_subjects_three_group(manifest, opt)

cat(paste0("\nTotal included subjects: ", length(uniq_sub)),
    file = logname, append = TRUE, sep = "\n")

# Intersect virus columns across all plates
common_virs <- Reduce(intersect, lapply(posfiles, colnames))
for (i in seq_along(plates))
    posfiles[[i]] <- posfiles[[i]][, common_virs]
common_virs <- common_virs[-c(1:3)]   # drop id, type, species header columns

cat(paste0("\nTotal included viruses: ", length(common_virs)),
    file = logname, append = TRUE, sep = "\n")

# ---------------------------------------------------------------------------
# Convert VirScan calls to binary (>1 hit threshold, same as original script)
# ---------------------------------------------------------------------------

cat("\nStart converting VirScan calls to binary calls",
    file = logname, append = TRUE, sep = "\n")

viral_calls           <- data.frame(matrix(0L, nrow = length(uniq_sub),
                                            ncol = length(common_virs) + 1))
colnames(viral_calls) <- c("ID", common_virs)
viral_calls$ID        <- uniq_sub

for (i in seq_len(nrow(viral_calls))) {
    print(paste("Looping:", i, ":", nrow(viral_calls)))

    id    <- viral_calls$ID[i]
    mp    <- manifest$plate[which(manifest$subject == id)][[1]]
    rloc  <- which(posfiles[[mp]][, 1] == id)[1]

    myCounts      <- data.frame(t(posfiles[[mp]][rloc,
                          which(colnames(posfiles[[mp]]) %in% common_virs)]))
    myCounts[, 1] <- as.numeric(myCounts[, 1])

    if (length(which((common_virs == row.names(myCounts)) == FALSE)) > 0) {
        cat(paste("\nError:", id, ": plate", mp,
                  ": Viruses out of order. This should never appear."),
            file = logname, append = TRUE, sep = "\n")
    }

    viral_calls[i, -1] <- ifelse(myCounts > 1, 1L, 0L)
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(posfiles)

# ---------------------------------------------------------------------------
# Build covariate data frame and model formula
# ---------------------------------------------------------------------------

datfile       <- build_datfile_three_group(uniq_sub, manifest, opt)
datfile$virus <- NA

cat(capture.output(print(head(datfile))), file = logname, append = TRUE, sep = "\n")

formula <- build_formula_three_group(datfile,
                                     predictor    = "virus",
                                     ignore_plate = opt$ignore_plate)
cat(paste("formula", formula), file = logname, append = TRUE, sep = "\n")

# ---------------------------------------------------------------------------
# Per-group counts for frequency calculations
# ---------------------------------------------------------------------------

n_g1 <- sum(datfile$group == opt$group1)
n_g2 <- sum(datfile$group == opt$group2)
n_g3 <- sum(datfile$group == opt$group3)

cat(paste0("\nTotal ", opt$group1, ": ", n_g1), file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal ", opt$group2, ": ", n_g2), file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal ", opt$group3, ": ", n_g3), file = logname, append = TRUE, sep = "\n")

# ---------------------------------------------------------------------------
# Initialise results table
# ---------------------------------------------------------------------------

g2_vs_g1 <- paste0(opt$group2, "_vs_", opt$group1)
g3_vs_g1 <- paste0(opt$group3, "_vs_", opt$group1)

result_cols <- c(
    "species",
    paste0("freq_", opt$group1),
    paste0("freq_", opt$group2),
    paste0("freq_", opt$group3),
    paste0("beta_", g2_vs_g1), paste0("se_", g2_vs_g1), paste0("pval_", g2_vs_g1),
    paste0("beta_", g3_vs_g1), paste0("se_", g3_vs_g1), paste0("pval_", g3_vs_g1)
)

pvalues           <- data.frame(matrix(NA, nrow = length(common_virs),
                                        ncol = length(result_cols)))
colnames(pvalues) <- result_cols
pvalues$species   <- common_virs

# ---------------------------------------------------------------------------
# Main loop: fit multinomial model per virus
# ---------------------------------------------------------------------------

cat("\nStart loop over virus multinomial regression analysis:",
    file = logname, append = TRUE, sep = "\n")

for (i in seq_along(common_virs)) {
    print(paste("Looping:", i, ":", length(common_virs)))

    datfile$virus <- viral_calls[, which(colnames(viral_calls) == common_virs[i])]

    n_g1_pos <- sum(datfile$virus == 1 & datfile$group == opt$group1)
    n_g2_pos <- sum(datfile$virus == 1 & datfile$group == opt$group2)
    n_g3_pos <- sum(datfile$virus == 1 & datfile$group == opt$group3)

    pvalues[i, paste0("freq_", opt$group1)] <- n_g1_pos / n_g1
    pvalues[i, paste0("freq_", opt$group2)] <- n_g2_pos / n_g2
    pvalues[i, paste0("freq_", opt$group3)] <- n_g3_pos / n_g3

    total_pos <- n_g1_pos + n_g2_pos + n_g3_pos
    total_n   <- n_g1 + n_g2 + n_g3

    if (total_pos %in% c(0, total_n)) {
        next
    }

    results <- multinom_reg(datfile, formula, "virus", logname, opt)

    pvalues[i, paste0("beta_", g2_vs_g1)] <- results$g2[1]
    pvalues[i, paste0("se_",   g2_vs_g1)] <- results$g2[2]
    pvalues[i, paste0("pval_", g2_vs_g1)] <- results$g2[3]

    pvalues[i, paste0("beta_", g3_vs_g1)] <- results$g3[1]
    pvalues[i, paste0("se_",   g3_vs_g1)] <- results$g3[2]
    pvalues[i, paste0("pval_", g3_vs_g1)] <- results$g3[3]
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

write_results_three_group(pvalues, output_base)

cat("\nAnalysis complete.", file = logname, append = TRUE, sep = "\n")
