#!/usr/bin/env Rscript

#   Three-group multinomial logistic regression for PhIP-seq peptide-level hits.
#
#   Extends Multi_Plate_Case_Control_Peptide_Regression.R from two groups to three
#   by replacing binary logistic regression (glm / binomial) with multinomial
#   logistic regression (nnet::multinom).  The original two-group scripts are
#   unchanged; this script lives alongside them.
#
#   Model
#   -----
#   For every peptide, fit:
#       group ~ peptide [+ age] [+ sex] [+ plate]
#   where group is a three-level factor with group1 as the reference.
#   This yields two sets of coefficients per peptide:
#       group2 vs. group1  (e.g. control1 vs. case)
#       group3 vs. group1  (e.g. control2 vs. case)
#
#   Output columns (per peptide, sorted by ascending min p-value):
#       peptide, species,
#       freq_<group1>, freq_<group2>, freq_<group3>,
#       beta_<group2>_vs_<group1>, se_<group2>_vs_<group1>, pval_<group2>_vs_<group1>,
#       beta_<group3>_vs_<group1>, se_<group3>_vs_<group1>, pval_<group3>_vs_<group1>
#
#   Usage example (parallel to the two-group script):
#       module load r
#       Multi_Plate_Three_Group_Peptide_Regression.R \
#           -z 5 --type "pemphigus serum" \
#           -a "PF Patient" -b "Non Endemic Control" -c "Healthy Control" \
#           --zfile_basename Zscores.csv \
#           -o /path/to/output \
#           -p /path/to/out.plate13 -p /path/to/out.plate14 \
#           -p /path/to/out.plate21 --sex F
#
#   Dependencies: argparse, data.table, nnet, fs
#   Sourced libs:  GenoLib.R  (read_zfile, read_multiple_manifests)
#                  ThreeGroupLib.R (build_datfile_three_group, select_subjects_three_group,
#                                   multinom_reg, build_formula_three_group,
#                                   init_log_three_group, write_results_three_group)

options(warn = 1)

library("argparse")
args       <- commandArgs()
scriptname <- sub("--file=", "", args[grepl("--file=", args)])
scriptdir  <- dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(paste(scriptdir, "GenoLib.R",        sep = "/"))
source(paste(scriptdir, "ThreeGroupLib.R",  sep = "/"))

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = scriptname)
parser$add_argument("--study",  type = "character", action = "append",
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
    help = "Z-score threshold for peptide positivity [default %(default)s]",
    metavar = "double")
parser$add_argument("-s", "--sex", type = "character", default = "",
    help = "limit sex", metavar = "sex")
parser$add_argument("-o", "--output_dir", type = "character", default = "./",
    help = "output directory [default %(default)s]", metavar = "directory")
parser$add_argument("-p", "--plate", type = "character", required = TRUE,
    action = "append",
    help = "plate directory (repeat for each plate)", metavar = "directory")
parser$add_argument("--zfile_basename", type = "character",
    default = "Zscores.csv",
    help = "Z-score file basename [default %(default)s]",
    metavar = "Zscores file basename")
parser$add_argument("--counts", action = "store_true",
    help = "values are raw counts, not Z-scores")
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
    paste("Multiplate_Three_Group_Peptide_Comparison",
          fs::path_ext_remove(basename(opt$zfile_basename)),
          paste(opt$study, collapse = "-"), opt$type,
          paste(groups_to_compare, collapse = "-"),
          "Z", Z, "sex", opt$sex, sep = "-")))

print(output_base)

logname <- paste0(output_base, ".log")

init_log_three_group(
    logname   = logname,
    analysis  = paste0(
        "Multi-plate three-group multinomial logistic regression for peptide ",
        "presence, adjusting for age, sex, and plate.\n",
        "Reference group (group1): ", opt$group1, "\n",
        "Z-score threshold = ", Z),
    plates    = plates,
    groups    = groups_to_compare
)

# ---------------------------------------------------------------------------
# Read Z-score files and manifests
# ---------------------------------------------------------------------------

Zfiles     <- list()
species_ids <- list()

for (i in seq_along(plates)) {
    results        <- read_zfile(paste(plates[i], opt$zfile_basename, sep = "/"))
    species_ids[[i]] <- results$species
    Zfiles[[i]]    <- results$zfile
}
rm(results)

manifest <- read_multiple_manifests(plates)
uniq_sub <- select_subjects_three_group(manifest, opt)

cat(paste0("\nUnique subjects: ",      paste(uniq_sub, collapse = ",")),
    file = logname, append = TRUE, sep = "\n")
cat(paste0("\nTotal included subjects: ", length(uniq_sub)),
    file = logname, append = TRUE, sep = "\n")

# Intersect peptides across all plates
common_peps <- Reduce(intersect, sapply(species_ids, `[`, 1))
print(paste("Total common peptides:", length(common_peps)))
cat(paste0("\nTotal included peptides: ", length(common_peps)),
    file = logname, append = TRUE, sep = "\n")

# Reorder / cull each Z-score matrix to common_peps
for (i in seq_along(plates)) {
    Zfiles[[i]] <- Zfiles[[i]][, c("id", common_peps)]
}

# ---------------------------------------------------------------------------
# Convert Z-scores to binary peptide calls (or keep as counts)
# ---------------------------------------------------------------------------

cat("\nStart converting Z scores to peptide binary calls",
    file = logname, append = TRUE, sep = "\n")

peptide_calls           <- data.frame(matrix(0L, nrow = length(uniq_sub),
                                              ncol = length(common_peps) + 1))
colnames(peptide_calls) <- c("ID", common_peps)
peptide_calls$ID        <- uniq_sub

for (i in seq_len(nrow(peptide_calls))) {
    print(paste("Looping:", i, ":", nrow(peptide_calls)))

    id   <- peptide_calls$ID[i]
    mp   <- manifest$plate[which(manifest$subject == id)][[1]]
    rlocs <- grep(id, Zfiles[[mp]][, 1])

    myZs <- data.frame(t(Zfiles[[mp]][c(1, rlocs),
                           which(Zfiles[[mp]][1, ] %in% common_peps)]))
    myZs[, c(2:3)] <- sapply(myZs[, c(2:3)], as.numeric)

    if (length(which((common_peps == myZs[, 1]) == FALSE)) > 0) {
        cat(paste("\nError:", id, ": plate", mp,
                  ": Peptides out of order. This should never appear."),
            file = logname, append = TRUE, sep = "\n")
    }

    mymins <- apply(myZs[, c(2:3)], 1, min)
    mymins[is.na(mymins)] <- 0

    if (!opt$counts) {
        peptide_calls[i, -1] <- ifelse(mymins > Z, 1L, 0L)
    } else {
        peptide_calls[i, -1] <- mymins
    }
}
cat("...Complete.", file = logname, append = TRUE, sep = "\n")

rm(Zfiles)

# ---------------------------------------------------------------------------
# Build covariate data frame and model formula
# ---------------------------------------------------------------------------

datfile         <- build_datfile_three_group(uniq_sub, manifest, opt)
datfile$peptide <- NA

cat(capture.output(print(head(datfile))), file = logname, append = TRUE, sep = "\n")

formula <- build_formula_three_group(datfile,
                                     predictor    = "peptide",
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
    "peptide", "species",
    paste0("freq_", opt$group1),
    paste0("freq_", opt$group2),
    paste0("freq_", opt$group3),
    paste0("beta_", g2_vs_g1), paste0("se_", g2_vs_g1), paste0("pval_", g2_vs_g1),
    paste0("beta_", g3_vs_g1), paste0("se_", g3_vs_g1), paste0("pval_", g3_vs_g1)
)

pvalues         <- data.frame(matrix(NA, nrow = length(common_peps),
                                      ncol = length(result_cols)))
colnames(pvalues) <- result_cols
pvalues$peptide   <- common_peps

# ---------------------------------------------------------------------------
# Main loop: fit multinomial model per peptide
# ---------------------------------------------------------------------------

cat("\nStart loop over peptide multinomial regression analysis:",
    file = logname, append = TRUE, sep = "\n")

for (i in seq_along(common_peps)) {
    print(paste("Looping:", i, ":", length(common_peps)))

    pvalues$species[i] <- species_ids[[1]]$species[
        which(species_ids[[1]]$id == common_peps[i])][1]

    datfile$peptide <- peptide_calls[, which(colnames(peptide_calls) == common_peps[i])]

    if (!opt$counts) {

        n_g1_pos <- sum(datfile$peptide == 1 & datfile$group == opt$group1)
        n_g2_pos <- sum(datfile$peptide == 1 & datfile$group == opt$group2)
        n_g3_pos <- sum(datfile$peptide == 1 & datfile$group == opt$group3)

        pvalues[i, paste0("freq_", opt$group1)] <- n_g1_pos / n_g1
        pvalues[i, paste0("freq_", opt$group2)] <- n_g2_pos / n_g2
        pvalues[i, paste0("freq_", opt$group3)] <- n_g3_pos / n_g3

        total_pos <- n_g1_pos + n_g2_pos + n_g3_pos
        total_n   <- n_g1 + n_g2 + n_g3

        # Skip peptides with no variation across all subjects
        if (total_pos %in% c(0, total_n)) {
            # Leave beta / se / pval as NA
            next
        }

    } else {
        # Log-transform raw counts (shift zeros to 0.001)
        datfile$peptide <- ifelse(datfile$peptide <= 0, 0.001, datfile$peptide)
        datfile$peptide <- log(datfile$peptide)

        pvalues[i, paste0("freq_", opt$group1)] <- "UNK"
        pvalues[i, paste0("freq_", opt$group2)] <- "UNK"
        pvalues[i, paste0("freq_", opt$group3)] <- "UNK"
    }

    results <- multinom_reg(datfile, formula, "peptide", logname, opt)

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
