
# ThreeGroupLib.R
#
# Shared helper functions for three-group multinomial PhIP-seq regression analyses.
# Mirrors the structure of GenoLib.R but replaces the binary logistic framework
# with multinomial logistic regression (nnet::multinom) and extends all helpers
# to accommodate a third group.
#
# Sourced by:
#   Multi_Plate_Three_Group_Peptide_Regression.R
#   Multi_Plate_Three_Group_VirHitFrac_Seropositivity_Regression.R
#   Multi_Plate_Three_Group_VirScan_Seropositivity_Regression.R
#
# Dependencies (in addition to those already loaded by the calling script):
#   data.table, nnet, fs

library(data.table)
library(nnet)


# ---------------------------------------------------------------------------
# build_datfile_three_group
#
# Builds the per-subject covariate data frame with a three-level factor
# for group membership.  Group1 is used as the reference level so that
# the two sets of multinomial coefficients represent:
#   group2 vs. group1  (control1 vs. case)
#   group3 vs. group1  (control2 vs. case)
#
# Arguments:
#   uniq_sub  : character vector of unique subject IDs (from select_subjects_three_group)
#   manifest  : aggregated manifest data frame (from read_multiple_manifests in GenoLib.R)
#   opt       : parsed argument list; must contain $group1, $group2, $group3
#
# Returns a data frame with columns: ID, group, sex, age, plate, lane
# ---------------------------------------------------------------------------

build_datfile_three_group <- function(uniq_sub, manifest, opt) {

    print("Building three-group datfile")

    # Lane is no longer used in the regression formula (it is nested within
    # plate, so the two are collinear).  We still create the column for
    # compatibility, but only populate it if the manifest actually has a
    # `lane` column -- otherwise it stays NA.
    has_lane <- "lane" %in% colnames(manifest)

    datfile <- data.frame(matrix(NA, nrow = length(uniq_sub),
                                 ncol = 6,
                                 dimnames = list(NULL, c("ID", "group", "sex", "age", "plate", "lane"))))
    datfile$ID <- uniq_sub

    for (i in seq_along(uniq_sub)) {
        print(paste("Looping:", i, ":", nrow(datfile)))
        man_loc <- which(manifest$subject == datfile$ID[i])[1]
        datfile$group[i] <- manifest$group[man_loc]
        datfile$age[i]   <- manifest$age[man_loc]
        datfile$sex[i]   <- manifest$sex[man_loc]
        datfile$plate[i] <- manifest$plate[man_loc]
        if (has_lane) {
            val <- manifest$lane[man_loc]
            datfile$lane[i] <- if (length(val) == 0) NA else val
        }
    }

    datfile$age   <- as.numeric(datfile$age)
    datfile$sex   <- as.factor(datfile$sex)
    datfile$plate <- as.factor(datfile$plate)
    datfile$lane  <- as.factor(datfile$lane)

    # Group is a factor; group1 is the reference level so coefficients are
    # interpreted as group2-vs-group1 and group3-vs-group1.
    datfile$group <- factor(datfile$group,
                            levels = c(opt$group1, opt$group2, opt$group3))

    print(head(datfile))
    return(datfile)
}


# ---------------------------------------------------------------------------
# select_subjects_three_group
#
# Extends GenoLib.R::select_subjects to filter for three groups and
# still apply the same optional sex / study / type filters.
#
# Arguments:
#   manifest : aggregated manifest data frame
#   opt      : parsed argument list; must contain $group1, $group2, $group3
#              and optionally $sex, $study, $type
#
# Returns a character vector of unique subject IDs that pass all filters.
# Quits (status 0) with a message if any group is empty after filtering.
# ---------------------------------------------------------------------------

select_subjects_three_group <- function(manifest, opt) {

    print("Selecting subjects (three-group)")

    all_groups <- c(opt$group1, opt$group2, opt$group3)

    uniq_sub <- unique(manifest$subject[which(manifest$group %in% all_groups)])

    if (opt$sex == "") {
        print("Sex is not set so not filtering on sex.")
    } else {
        print(paste0("Sex is set to ", opt$sex, ". Filtering."))
        uniq_sub <- intersect(uniq_sub,
                              unique(manifest$subject[which(manifest$sex == opt$sex)]))
    }

    if (length(opt$study) == 0) {
        print("Study is not set so not filtering on study.")
    } else {
        print(paste0("Study is set to ", paste(opt$study, collapse = " or "), ". Filtering."))
        uniq_sub <- intersect(uniq_sub,
                              unique(manifest$subject[which(manifest$study %in% opt$study)]))
    }

    if (opt$type == "") {
        print("Type is not set so not filtering on type.")
    } else {
        print(paste0("Type is set to ", opt$type, ". Filtering."))
        uniq_sub <- intersect(uniq_sub,
                              unique(manifest$subject[which(manifest$type == opt$type)]))
    }

    print(paste("Total subjects after filtering:", length(uniq_sub)))

    # Verify each group is non-empty
    for (grp in all_groups) {
        grp_subs <- intersect(unique(manifest$subject[which(manifest$group == grp)]), uniq_sub)
        print(paste0("  ", grp, ": n=", length(grp_subs)))
        if (length(grp_subs) == 0) {
            print(paste("Group", grp, "has no subjects after filtering. Quitting."))
            quit(save = "no", status = 0, runLast = FALSE)
        }
    }

    return(uniq_sub)
}


# ---------------------------------------------------------------------------
# multinom_reg
#
# Fits a multinomial logistic regression model and extracts the coefficients
# for the predictor of interest for each contrast (group2 vs. group1 and
# group3 vs. group1).
#
# The model formula is passed in as a string, e.g.
#   "group ~ predictor + age + sex + plate"
# where the response is always 'group' (a three-level factor with group1
# as the reference level set in build_datfile_three_group).
#
# Arguments:
#   df        : data frame from build_datfile_three_group plus the predictor column
#   formula   : model formula string
#   label     : name of the predictor column in df (e.g. "peptide" or "virus")
#   logname   : path to the log file for appending model summaries
#   opt       : parsed argument list (for group labels in column naming)
#
# Returns a named list with two elements, one per non-reference group:
#   $g2  : c(beta, se, pval) for group2 vs. group1
#   $g3  : c(beta, se, pval) for group3 vs. group1
# Each element is c(NA, NA, NA) if the label is absent from the model output
# (e.g. due to singularities).
# ---------------------------------------------------------------------------

multinom_reg <- function(df, formula, label, logname, opt) {

    # nnet::multinom prints iteration traces; suppress them in the log
    suppressMessages(
        fit <- nnet::multinom(as.formula(formula), data = df, trace = FALSE)
    )

    go <- summary(fit)
    cat(capture.output(print(go)), file = logname, append = TRUE, sep = "\n")

    # z-scores and two-sided p-values (multinom does not compute these natively)
    z_scores <- go$coefficients / go$standard.errors
    p_matrix  <- 2 * (1 - pnorm(abs(z_scores)))

    # Row names in multinom output are the non-reference group labels
    g2_label <- opt$group2
    g3_label <- opt$group3

    extract_row <- function(group_label) {
        if (!(group_label %in% rownames(go$coefficients)) ||
            !(label %in% colnames(go$coefficients))) {
            return(c(NA_real_, NA_real_, NA_real_))
        }
        beta <- go$coefficients[group_label, label]
        se   <- go$standard.errors[group_label, label]
        pval <- p_matrix[group_label, label]
        c(beta, se, pval)
    }

    return(list(g2 = extract_row(g2_label),
                g3 = extract_row(g3_label)))
}


# ---------------------------------------------------------------------------
# build_formula_three_group
#
# Constructs the multinomial model formula string, adding covariates only
# when they vary across the analysis sample (mirrors the logic in the
# original two-group scripts).
#
# Arguments:
#   datfile     : data frame from build_datfile_three_group
#   predictor   : name of the predictor column, e.g. "peptide" or "virus"
#   ignore_plate: logical; if TRUE, plate is excluded from the formula
#
# Returns a character string, e.g. "group ~ peptide + age + sex + plate"
# ---------------------------------------------------------------------------

build_formula_three_group <- function(datfile, predictor, ignore_plate = FALSE) {

    formula <- paste("group ~", predictor)

    if (length(unique(datfile$age)) > 1)
        formula <- paste(formula, "age", sep = " + ")
    if (length(unique(datfile$sex)) > 1)
        formula <- paste(formula, "sex", sep = " + ")
    if (length(unique(datfile$plate)) > 1 && !ignore_plate)
        formula <- paste(formula, "plate", sep = " + ")
    # Lane is not included: each lane contains 2 plates, so lane is nested
    # within (and largely redundant with) the plate covariate.  Including both
    # would risk rank-deficiency / collinearity issues.

    print(paste("formula", formula))
    return(formula)
}


# ---------------------------------------------------------------------------
# init_log_three_group
#
# Opens the log file and writes the standard header block (analysis type,
# z-score threshold, plates, groups).
#
# Arguments:
#   logname     : full path to the log file
#   analysis    : short description string, printed at the top of the log
#   plates      : character vector of plate directories
#   groups      : character vector of the three group labels
#   extra_lines : optional character vector of additional lines to append
#                 after the groups block (e.g. "Z-score threshold = 5")
# ---------------------------------------------------------------------------

init_log_three_group <- function(logname, analysis, plates, groups,
                                 extra_lines = character(0)) {

    cat(analysis, file = logname, sep = "\n")

    if (length(extra_lines) > 0)
        cat(extra_lines, file = logname, append = TRUE, sep = "\n")

    cat("\nPlates used in this analysis:", file = logname, append = TRUE, sep = "\n")
    for (p in plates)
        cat(p, file = logname, append = TRUE, sep = "\n")
    cat("\n", file = logname, append = TRUE)

    cat("\nGroups compared in this analysis:", file = logname, append = TRUE, sep = "\n")
    for (g in groups)
        cat(g, file = logname, append = TRUE, sep = "\n")
    cat("\n", file = logname, append = TRUE)
}


# ---------------------------------------------------------------------------
# write_results_three_group
#
# Writes the final results data frame to a CSV, sorted by ascending
# omnibus p-value (smallest of the two contrast p-values per feature),
# then by the group2-vs-group1 p-value as a tiebreaker.
#
# Arguments:
#   pvalues     : results data frame (columns defined by calling script)
#   output_base : path prefix; ".csv" is appended automatically
# ---------------------------------------------------------------------------

write_results_three_group <- function(pvalues, output_base) {

    # Determine the omnibus sort key: min p-value across the two contrasts
    pval_cols <- grep("^pval_", colnames(pvalues), value = TRUE)

    if (length(pval_cols) >= 2) {
        pvalues$pval_min <- apply(pvalues[, pval_cols], 1,
                                  function(x) min(x, na.rm = TRUE))
        pvalues$pval_min[!is.finite(pvalues$pval_min)] <- NA
        sort_col <- "pval_min"
    } else {
        sort_col <- pval_cols[1]
    }

    out <- pvalues[order(pvalues[[sort_col]], decreasing = FALSE, na.last = TRUE), ]

    # Drop the helper sort column before writing
    if ("pval_min" %in% colnames(out))
        out$pval_min <- NULL

    write.table(out, paste0(output_base, ".csv"),
                col.names = TRUE, sep = ",", row.names = FALSE, quote = FALSE)
}


print("Loaded ThreeGroupLib")
