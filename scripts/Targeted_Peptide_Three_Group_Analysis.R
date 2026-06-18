#!/usr/bin/env Rscript

#   Targeted_Peptide_Three_Group_Analysis.R
#
#   Stand-alone exploratory analysis for a user-supplied subset of peptide IDs
#   (e.g. tiles with homology to Dsg1 / Dsg3 or other antigens of interest).
#
#   Produces, in one HTML report:
#       1. Per-peptide multinomial regression restricted to the subset
#       2. Per-subject composite scores:
#            (a) count of subset peptides hit (above Z threshold)
#            (b) mean Z-score across subset
#          with Kruskal-Wallis and pairwise Wilcoxon tests across the 3 groups
#       3. Per-subject heatmap of Z-scores across the subset (rows = subjects
#          ordered by group, columns = peptides)
#
#   The script is intentionally independent of the genome-wide pipeline -- it
#   does not consume regression CSVs, but reads the raw Z-score files directly
#   so the composite scores are computed faithfully (same min-of-duplicates
#   per-plate logic).  This keeps the analysis honest and self-contained.
#
#   Usage:
#       module load r
#       Targeted_Peptide_Three_Group_Analysis.R \
#           -a "PF Patient" \
#           -b "Endemic Control" \
#           -c "Non Endemic Control" \
#           --antigen Dsg1 \
#           --peptide_list /path/to/dsg1_tile_ids.txt \
#           -z 5 \
#           -p out.plate13 -p out.plate14 [ -p ... ] \
#           --zfile_basename Zscores.csv \
#           -o /path/to/output_dir \
#           [--sex F] [--type "pemphigus serum"] [--study STUDY] \
#           [--ignore_plate]
#
#   --peptide_list is a plain text file with one tile ID per line (header
#   optional; non-numeric lines are skipped silently).
#
#   Output files (named by antigen, Z, sex):
#       Targeted_<antigen>_<groups>_Z<z>_sex<s>.csv     -- per-peptide results
#       Targeted_<antigen>_<groups>_Z<z>_sex<s>_composite.csv -- per-subject scores
#       Targeted_<antigen>_<groups>_Z<z>_sex<s>.html    -- standalone HTML report
#
#   Dependencies: argparse, data.table, nnet, fs, ggplot2, dplyr, tidyr,
#                 rmarkdown, knitr, kableExtra, scales
#   Sourced libs: GenoLib.R       (read_zfile, read_multiple_manifests)
#                 ThreeGroupLib.R (build_datfile_three_group,
#                                  select_subjects_three_group, multinom_reg,
#                                  build_formula_three_group)

options(warn = 1)

library("argparse")
args      <- commandArgs()
scriptdir <- dirname(sub("--file=", "", args[grepl("--file=", args)]))
source(file.path(scriptdir, "GenoLib.R"))
source(file.path(scriptdir, "ThreeGroupLib.R"))

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parser <- ArgumentParser()
parser$add_argument("-a", "--group1", required = TRUE,
    help = "first group label  (reference / cases)", metavar = "group")
parser$add_argument("-b", "--group2", required = TRUE,
    help = "second group label (control 1)", metavar = "group")
parser$add_argument("-c", "--group3", required = TRUE,
    help = "third group label  (control 2)", metavar = "group")
parser$add_argument("--antigen", required = TRUE,
    help = "antigen name for output labelling (e.g. Dsg1)", metavar = "name")
parser$add_argument("--peptide_list", required = TRUE,
    help = "path to text file containing one tile ID per line",
    metavar = "file")
parser$add_argument("-z", "--zscore", type = "double", default = 5,
    help = "Z-score threshold for hit calling [default %(default)s]",
    metavar = "double")
parser$add_argument("-p", "--plate", required = TRUE, action = "append",
    help = "plate directory (repeat for each plate)", metavar = "directory")
parser$add_argument("--zfile_basename", default = "Zscores.csv",
    help = "Z-score file basename [default %(default)s]", metavar = "name")
parser$add_argument("--study", action = "append",
    help = "limit study", metavar = "study")
parser$add_argument("-t", "--type", default = "",
    help = "limit type", metavar = "type")
parser$add_argument("-s", "--sex", default = "",
    help = "limit sex", metavar = "sex")
parser$add_argument("-o", "--output_dir", default = "./",
    help = "output directory [default %(default)s]", metavar = "directory")
parser$add_argument("--ignore_plate", action = "store_true",
    help = "exclude plate from the regression formula")
opt <- parser$parse_args()

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Output naming
# ---------------------------------------------------------------------------

groups_to_compare <- c(opt$group1, opt$group2, opt$group3)
sex_tag <- if (nchar(opt$sex) == 0) "both" else opt$sex

output_base <- file.path(opt$output_dir,
    paste0("Targeted_", gsub(" ", "_", opt$antigen), "_",
           gsub(" ", "_", paste(groups_to_compare, collapse = "-")),
           "_Z", opt$zscore, "_sex", sex_tag))

logname <- paste0(output_base, ".log")
cat("Targeted peptide analysis\n", file = logname)
cat(paste("Antigen:", opt$antigen, "\n"), file = logname, append = TRUE)
cat(paste("Z-score threshold:", opt$zscore, "\n"), file = logname, append = TRUE)
cat(paste("Sex filter:", ifelse(nchar(opt$sex) == 0, "(both)", opt$sex), "\n"),
    file = logname, append = TRUE)
cat(paste("Reference group:", opt$group1, "\n"), file = logname, append = TRUE)
cat(paste("Plates:", paste(opt$plate, collapse = ", "), "\n"),
    file = logname, append = TRUE)

# ---------------------------------------------------------------------------
# Read the peptide ID list (one ID per line; skip non-numeric tokens)
# ---------------------------------------------------------------------------

raw <- readLines(opt$peptide_list)
raw <- trimws(raw)
raw <- raw[nchar(raw) > 0]
# Drop any header-like row (non-numeric)
target_ids <- raw[grepl("^[0-9]+$", raw)]
n_dropped  <- length(raw) - length(target_ids)
if (n_dropped > 0)
  cat(paste("Skipped", n_dropped, "non-numeric line(s) from peptide list\n"),
      file = logname, append = TRUE)
if (length(target_ids) == 0)
  stop("Peptide list is empty after parsing. Check the file format.")

cat(paste("Target peptides loaded:", length(target_ids), "\n"),
    file = logname, append = TRUE)

# ---------------------------------------------------------------------------
# Load Z-score data and manifests across plates
# ---------------------------------------------------------------------------

plates <- opt$plate
Zfiles      <- list()
species_ids <- list()
for (i in seq_along(plates)) {
  results          <- read_zfile(file.path(plates[i], opt$zfile_basename))
  species_ids[[i]] <- results$species
  Zfiles[[i]]      <- results$zfile
}
rm(results)

manifest <- read_multiple_manifests(plates)
uniq_sub <- select_subjects_three_group(manifest, opt)

# Restrict target IDs to those present on ALL plates
common_peps <- Reduce(intersect, sapply(species_ids, `[`, 1))
present_ids <- intersect(target_ids, common_peps)
missing_ids <- setdiff(target_ids, common_peps)
cat(paste("Target peptides present on all plates:", length(present_ids), "\n"),
    file = logname, append = TRUE)
if (length(missing_ids) > 0) {
  cat(paste("Target peptides MISSING from at least one plate:\n  ",
            paste(missing_ids, collapse = ", "), "\n"),
      file = logname, append = TRUE)
}
if (length(present_ids) == 0)
  stop("None of the target peptide IDs were found in the Z-score files.")

# Map IDs to species labels for downstream display
sp_lookup <- species_ids[[1]]
sp_for_id <- function(id) {
  hit <- sp_lookup$species[sp_lookup$id == id]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# ---------------------------------------------------------------------------
# Build per-subject Z-score and hit-call matrices for the target peptides only
# ---------------------------------------------------------------------------

Z <- opt$zscore

zmat   <- matrix(NA_real_, nrow = length(uniq_sub),
                 ncol = length(present_ids),
                 dimnames = list(uniq_sub, present_ids))
hitmat <- matrix(0L, nrow = length(uniq_sub),
                 ncol = length(present_ids),
                 dimnames = list(uniq_sub, present_ids))

for (i in seq_along(uniq_sub)) {
  id   <- uniq_sub[i]
  mp   <- manifest$plate[which(manifest$subject == id)][[1]]
  rlocs <- grep(id, Zfiles[[mp]][, 1])

  # Same min-of-duplicates extraction as the genome-wide script
  cols_to_keep <- which(Zfiles[[mp]][1, ] %in% present_ids)
  myZs <- data.frame(t(Zfiles[[mp]][c(1, rlocs), cols_to_keep]))
  myZs[, -1] <- sapply(myZs[, -1, drop = FALSE], as.numeric)

  # Reorder to match present_ids exactly
  ord <- match(present_ids, myZs[, 1])
  ordered_vals <- myZs[ord, -1, drop = FALSE]

  if (ncol(ordered_vals) == 1) {
    # Only one duplicate row available
    z_per_peptide <- as.numeric(ordered_vals[, 1])
  } else {
    z_per_peptide <- apply(ordered_vals, 1, function(v) {
      v <- v[!is.na(v)]
      if (length(v) == 0) NA_real_ else min(v)
    })
  }

  zmat[i, ]   <- z_per_peptide
  hitmat[i, ] <- ifelse(!is.na(z_per_peptide) & z_per_peptide > Z, 1L, 0L)
}

rm(Zfiles)

# ---------------------------------------------------------------------------
# Per-peptide multinomial regression on the target subset
# ---------------------------------------------------------------------------

datfile         <- build_datfile_three_group(uniq_sub, manifest, opt)
datfile$peptide <- NA_integer_

formula <- build_formula_three_group(datfile, predictor = "peptide",
                                     ignore_plate = opt$ignore_plate)
cat(paste("Formula:", formula, "\n"), file = logname, append = TRUE)

g2_vs_g1 <- paste0(opt$group2, "_vs_", opt$group1)
g3_vs_g1 <- paste0(opt$group3, "_vs_", opt$group1)

per_pep <- data.frame(
  peptide                                  = present_ids,
  species                                  = vapply(present_ids, sp_for_id, character(1)),
  n_hit                                    = colSums(hitmat),
  freq                                     = colMeans(hitmat),
  stringsAsFactors                         = FALSE
)
per_pep[[paste0("freq_", opt$group1)]] <- NA_real_
per_pep[[paste0("freq_", opt$group2)]] <- NA_real_
per_pep[[paste0("freq_", opt$group3)]] <- NA_real_
per_pep[[paste0("beta_", g2_vs_g1)]]   <- NA_real_
per_pep[[paste0("se_",   g2_vs_g1)]]   <- NA_real_
per_pep[[paste0("pval_", g2_vs_g1)]]   <- NA_real_
per_pep[[paste0("beta_", g3_vs_g1)]]   <- NA_real_
per_pep[[paste0("se_",   g3_vs_g1)]]   <- NA_real_
per_pep[[paste0("pval_", g3_vs_g1)]]   <- NA_real_

n_g1 <- sum(datfile$group == opt$group1)
n_g2 <- sum(datfile$group == opt$group2)
n_g3 <- sum(datfile$group == opt$group3)

cat(paste("Per-peptide regression over", length(present_ids), "peptides\n"),
    file = logname, append = TRUE)

for (i in seq_along(present_ids)) {
  pep_id          <- present_ids[i]
  datfile$peptide <- hitmat[, pep_id]

  pos_g1 <- sum(datfile$peptide == 1 & datfile$group == opt$group1)
  pos_g2 <- sum(datfile$peptide == 1 & datfile$group == opt$group2)
  pos_g3 <- sum(datfile$peptide == 1 & datfile$group == opt$group3)
  per_pep[i, paste0("freq_", opt$group1)] <- pos_g1 / n_g1
  per_pep[i, paste0("freq_", opt$group2)] <- pos_g2 / n_g2
  per_pep[i, paste0("freq_", opt$group3)] <- pos_g3 / n_g3

  total_pos <- pos_g1 + pos_g2 + pos_g3
  if (total_pos %in% c(0, nrow(datfile))) next

  res <- multinom_reg(datfile, formula, "peptide", logname, opt)

  per_pep[i, paste0("beta_", g2_vs_g1)] <- res$g2[1]
  per_pep[i, paste0("se_",   g2_vs_g1)] <- res$g2[2]
  per_pep[i, paste0("pval_", g2_vs_g1)] <- res$g2[3]
  per_pep[i, paste0("beta_", g3_vs_g1)] <- res$g3[1]
  per_pep[i, paste0("se_",   g3_vs_g1)] <- res$g3[2]
  per_pep[i, paste0("pval_", g3_vs_g1)] <- res$g3[3]
}

# Apply BH FDR within the subset (small set, so this is conservative-friendly)
per_pep$pval_g2_fdr <- p.adjust(per_pep[[paste0("pval_", g2_vs_g1)]],
                                method = "BH")
per_pep$pval_g3_fdr <- p.adjust(per_pep[[paste0("pval_", g3_vs_g1)]],
                                method = "BH")
per_pep$pval_min    <- pmin(per_pep[[paste0("pval_", g2_vs_g1)]],
                            per_pep[[paste0("pval_", g3_vs_g1)]],
                            na.rm = TRUE)
per_pep <- per_pep[order(per_pep$pval_min, na.last = TRUE), ]

write.table(per_pep, paste0(output_base, ".csv"),
            row.names = FALSE, sep = ",", quote = FALSE)

# ---------------------------------------------------------------------------
# Per-subject composite scores
# ---------------------------------------------------------------------------

comp <- data.frame(
  subject   = uniq_sub,
  group     = as.character(datfile$group),
  n_hit     = rowSums(hitmat, na.rm = TRUE),
  mean_z    = rowMeans(zmat, na.rm = TRUE),
  max_z     = apply(zmat, 1, function(v) {
                  v <- v[!is.na(v)]; if (length(v) == 0) NA_real_ else max(v)
              }),
  stringsAsFactors = FALSE
)
comp$frac_hit <- comp$n_hit / length(present_ids)

write.table(comp, paste0(output_base, "_composite.csv"),
            row.names = FALSE, sep = ",", quote = FALSE)

# ---------------------------------------------------------------------------
# Save intermediate objects for the Rmd to consume
# ---------------------------------------------------------------------------

rds_path <- paste0(output_base, "_data.rds")
saveRDS(list(
  per_pep      = per_pep,
  comp         = comp,
  zmat         = zmat,
  hitmat       = hitmat,
  present_ids  = present_ids,
  missing_ids  = missing_ids,
  groups       = groups_to_compare,
  group1       = opt$group1,
  group2       = opt$group2,
  group3       = opt$group3,
  antigen      = opt$antigen,
  zscore       = Z,
  sex          = opt$sex,
  output_base  = output_base
), rds_path)

# ---------------------------------------------------------------------------
# Render the report (Rmd lives alongside this script)
# ---------------------------------------------------------------------------

rmd_path <- file.path(scriptdir, "Targeted_Peptide_Three_Group_Report.Rmd")
if (!file.exists(rmd_path))
  stop(paste("Rmd not found:", rmd_path))

if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("Package 'rmarkdown' is required.")

html_path <- paste0(output_base, ".html")
cat("\nRendering HTML report -> ", html_path, "\n")

rmarkdown::render(
  input       = rmd_path,
  output_file = html_path,
  params      = list(rds_path = rds_path),
  envir       = new.env(parent = globalenv()),
  quiet       = FALSE
)

cat("\nOutputs:\n")
cat("  CSV (per peptide):  ", paste0(output_base, ".csv"), "\n")
cat("  CSV (per subject):  ", paste0(output_base, "_composite.csv"), "\n")
cat("  RDS (raw matrices): ", rds_path, "\n")
cat("  HTML report:        ", html_path, "\n")
cat("  Log:                ", logname, "\n")
