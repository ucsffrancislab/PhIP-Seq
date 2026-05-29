#!/usr/bin/env Rscript

#   Three-group version of Count_Viral_Tile_Hit_Fraction.R
#
#   Part 1 — Viral_Frac_Hits
#   -------------------------
#   Creates a per-plate "Viral_Frac_Hits_*.csv" table containing the fraction
#   of tiles hit for each virus for each sample (identical logic to the original),
#   but now selects subjects from all three groups and encodes all three group
#   labels in the output filename so that
#   Multi_Plate_Three_Group_VirHitFrac_Seropositivity_Regression.R can locate it.
#
#   Part 2 — Viral_Sero (exploratory, single-plate)
#   -------------------------------------------------
#   Replaces the two-group prop.test with a chi-square test over a 3×2
#   contingency table (seropositive / seronegative × three groups).
#   Reports a single omnibus p-value per virus plus the per-group seropositive
#   fraction, ordered by ascending p-value.
#
#   Output filenames:
#       Viral_Frac_Hits_<zfile_stem>-<study>-<type>-<g1>-<g2>-<g3>-Z-<Z>-sex-<sex>.csv
#       Viral_Sero_<zfile_stem>-<study>-<type>-<g1>-<g2>-<g3>-Z-<Z>-sex-<sex>-Vir_hit_frac-<vf>.csv
#
#   Usage example (mirrors two-group script):
#       module load r
#       Count_Viral_Tile_Hit_Fraction_Three_Group.R \
#           -a "PF Patient" -b "Non Endemic Control" -c "Healthy Control" \
#           --type "pemphigus serum" \
#           -z 5 \
#           --manifest /path/to/out.plate13/manifest.csv \
#           --output_dir /path/to/out.plate13 \
#           --zfilename /path/to/out.plate13/Zscores.csv
#
#   Dependencies: argparse, data.table, fs
#   Sourced libs:  GenoLib.R       (read_zfile, select_subjects)
#                  ThreeGroupLib.R (select_subjects_three_group)

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
parser$add_argument("-s", "--sex", type = "character", default = "",
    help = "limit sex", metavar = "sex")
parser$add_argument("--study", type = "character", action = "append",
    help = "limit study", metavar = "study")
parser$add_argument("-t", "--type", type = "character", default = "",
    help = "limit type", metavar = "type")
parser$add_argument("-a", "--group1", type = "character", required = TRUE,
    help = "first group", metavar = "group")
parser$add_argument("-b", "--group2", type = "character", required = TRUE,
    help = "second group", metavar = "group")
parser$add_argument("-c", "--group3", type = "character", required = TRUE,
    help = "third group", metavar = "group")
parser$add_argument("-z", "--zscore", type = "double", default = 3.5,
    help = "Z-score threshold [default %(default)s]", metavar = "double")
parser$add_argument("-m", "--manifest", type = "character", default = NULL,
    required = TRUE,
    help = "manifest file path", metavar = "manifest")
parser$add_argument("--output_dir", type = "character", default = "./",
    help = "output directory [default %(default)s]", metavar = "directory")
parser$add_argument("--zfilename", type = "character", default = "out/Zscores.csv",
    help = "path to Zscores file [default %(default)s]", metavar = "Zscores file")
opt <- parser$parse_args()

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

groups_to_compare <- c(opt$group1, opt$group2, opt$group3)
print("Comparing these groups"); print(groups_to_compare)

Z <- opt$zscore

library(data.table)

# ---------------------------------------------------------------------------
# Read manifest, Z-score file, and select subjects across all three groups
# ---------------------------------------------------------------------------

print("Read in the manifest file")
manifest <- data.frame(data.table::fread(opt$manifest, sep = ",", header = TRUE))

results    <- read_zfile(opt$zfilename)
species_id <- results$species
Zfile      <- results$zfile
rm(results)

print("Unique samples to keep (all three groups)")
uniq_sub <- select_subjects_three_group(manifest, opt)

# Subset Z-score matrix to the selected subjects (keep original + dup rows)
to_keep <- 1
for (u in uniq_sub) {
    possible_ids <- grep(u, Zfile[, 1])
    mids  <- Zfile[possible_ids, 1]
    locs  <- grepl("dup", mids)
    to_keep <- c(to_keep,
                 possible_ids[c(which(!locs)[1], which(locs)[1])])
}
Zfile1 <- Zfile[to_keep, ]
rm(Zfile)

# ---------------------------------------------------------------------------
# Part 1: Compute viral tile hit fractions
# ---------------------------------------------------------------------------

print("Unique species (column name is 'species' NOT 'Species')")
uniq_spec <- unique(species_id$species)
print(uniq_spec[1:5])

print("Shell file for viral fractions")
virfracs           <- data.frame(matrix(0, nrow = length(uniq_sub),
                                         ncol = length(uniq_spec) + 1))
colnames(virfracs) <- c("id", uniq_spec)
virfracs$id        <- uniq_sub

print(paste("Looping over", ncol(virfracs), "columns"))
for (j in seq(2, ncol(virfracs))) {
    print(paste("Looping", j, ":", ncol(virfracs)))
    sp      <- colnames(virfracs)[j]
    vir_ids <- species_id$id[which(species_id$species == sp)]
    sfile   <- Zfile1[, c(1, which(Zfile1[1, ] %in% vir_ids))]

    for (i in seq_len(nrow(virfracs))) {
        vf <- sfile[grep(virfracs$id[i], sfile[, 1]), ]
        if (ncol(vf) > 2) {
            n_positive <- length(which(
                as.numeric(apply(vf[, -1], 2, min)) > Z))
        } else {
            n_positive <- ifelse(as.numeric(min(vf[, 2])) > Z, 1, 0)
        }
        virfracs[i, j] <- n_positive / length(vir_ids)
    }
}

# Output filename encodes all three groups so the multiplate script can
# locate this file by reconstructing the same name.
outfile_virfracs <- paste0(
    opt$output_dir, "/",
    gsub(" ", "_",
         paste("Viral_Frac_Hits",
               fs::path_ext_remove(basename(opt$zfilename)),
               paste(opt$study, collapse = "-"), opt$type,
               paste(groups_to_compare, collapse = "-"),
               "Z", Z,
               "sex", opt$sex,
               sep = "-")),
    ".csv")

print(paste0("Writing ", outfile_virfracs))
write.table(virfracs, outfile_virfracs,
            col.names = TRUE, sep = ",", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Part 2: Single-plate chi-square seropositivity comparison (three groups)
#
# For each virus, build a 3x2 contingency table:
#   rows  = group1, group2, group3
#   cols  = seropositive (hit fraction > Vir_frac), seronegative
# and run a chi-square test for independence.
# ---------------------------------------------------------------------------

Vir_frac <- 0.05

virfracs_sub <- virfracs[which(virfracs$id %in% uniq_sub), ]
print(paste("Subjects in virfracs:", nrow(virfracs_sub)))

# Per-group subject vectors
g1_subs <- intersect(unique(manifest$subject[manifest$group == opt$group1]), uniq_sub)
g2_subs <- intersect(unique(manifest$subject[manifest$group == opt$group2]), uniq_sub)
g3_subs <- intersect(unique(manifest$subject[manifest$group == opt$group3]), uniq_sub)

print(paste("Group1:", length(g1_subs),
            " Group2:", length(g2_subs),
            " Group3:", length(g3_subs)))

print("Create results table for three-group chi-square seropositivity test")
pvalues           <- data.frame(matrix(NA, nrow = ncol(virfracs) - 1, ncol = 5))
colnames(pvalues) <- c("species",
                       paste0("freq_", opt$group1),
                       paste0("freq_", opt$group2),
                       paste0("freq_", opt$group3),
                       "pval")
pvalues$species <- colnames(virfracs)[-1]

for (i in seq_len(nrow(pvalues))) {
    print(paste("Looping:", i, ":", nrow(pvalues)))
    sp      <- pvalues$species[i]
    sp_col  <- which(colnames(virfracs_sub) == sp)

    frac_g1 <- as.numeric(virfracs_sub[virfracs_sub$id %in% g1_subs, sp_col])
    frac_g2 <- as.numeric(virfracs_sub[virfracs_sub$id %in% g2_subs, sp_col])
    frac_g3 <- as.numeric(virfracs_sub[virfracs_sub$id %in% g3_subs, sp_col])

    pos_g1 <- sum(frac_g1 > Vir_frac, na.rm = TRUE)
    pos_g2 <- sum(frac_g2 > Vir_frac, na.rm = TRUE)
    pos_g3 <- sum(frac_g3 > Vir_frac, na.rm = TRUE)

    n_g1 <- length(g1_subs)
    n_g2 <- length(g2_subs)
    n_g3 <- length(g3_subs)

    pvalues[i, paste0("freq_", opt$group1)] <- pos_g1 / n_g1
    pvalues[i, paste0("freq_", opt$group2)] <- pos_g2 / n_g2
    pvalues[i, paste0("freq_", opt$group3)] <- pos_g3 / n_g3

    # 3x2 contingency table: rows = groups, cols = pos / neg
    ctable <- matrix(
        c(pos_g1, n_g1 - pos_g1,
          pos_g2, n_g2 - pos_g2,
          pos_g3, n_g3 - pos_g3),
        nrow = 3, byrow = TRUE,
        dimnames = list(
            c(opt$group1, opt$group2, opt$group3),
            c("seropositive", "seronegative")
        )
    )

    # Chi-square test requires at least one cell with an observed count;
    # if everyone is pos or neg across all groups, p-value is NA.
    total_pos <- pos_g1 + pos_g2 + pos_g3
    total_n   <- n_g1 + n_g2 + n_g3

    if (total_pos %in% c(0, total_n)) {
        pvalues$pval[i] <- NA
    } else {
        # suppress.warnings: chisq.test warns on cells with expected < 5,
        # which is common for rare viruses; the p-value is still returned.
        ct_result       <- suppressWarnings(chisq.test(ctable, correct = FALSE))
        pvalues$pval[i] <- ct_result$p.value
    }
}

outfile_sero <- paste0(
    opt$output_dir, "/",
    gsub(" ", "_",
         paste("Viral_Sero",
               fs::path_ext_remove(basename(opt$zfilename)),
               paste(opt$study, collapse = "-"), opt$type,
               paste(groups_to_compare, collapse = "-"),
               "Z", Z,
               "sex", opt$sex,
               "Vir_hit_frac", Vir_frac,
               sep = "-")),
    ".csv")

print(paste0("Writing ", outfile_sero))
write.table(
    pvalues[order(pvalues$pval, decreasing = FALSE, na.last = TRUE), ],
    outfile_sero,
    col.names = TRUE, sep = ",", row.names = FALSE, quote = FALSE)
