#!/usr/bin/env Rscript

#   Render_Three_Group_Combined_Report.R
#
#   Generates one combined cross-condition HTML report covering the full grid
#   of Z-score thresholds and sex strata.  By default it also renders each
#   individual per-condition report into a subdirectory, so the combined
#   report's "per-condition reports" section can link to them.
#
#   Usage:
#       module load r
#       Render_Three_Group_Combined_Report.R \
#           -a "PF Patient" \
#           -b "Endemic Control" \
#           -c "Non Endemic Control" \
#           -z "3.5,5,10" \
#           -s ",M,F" \
#           -o /path/to/output_dir \
#           [--fdr] [--pval 0.05] [--top_n 20] \
#           [--skip_per_condition] [--per_condition_subdir per_condition_reports]
#
#   The --zscores and --sexes flags take comma-separated lists.  For sex,
#   an empty token means "both sexes" — e.g. ",M,F" means: both, male only,
#   female only.
#
#   Dependencies: rmarkdown, knitr, plus all packages listed in the Rmds.

library("argparse")

args      <- commandArgs()
scriptdir <- dirname(sub("--file=", "", args[grepl("--file=", args)]))

parser <- ArgumentParser()
parser$add_argument("-a", "--group1", required = TRUE, metavar = "group",
    help = "first group label (reference / cases)")
parser$add_argument("-b", "--group2", required = TRUE, metavar = "group",
    help = "second group label (control 1)")
parser$add_argument("-c", "--group3", required = TRUE, metavar = "group",
    help = "third group label (control 2)")
parser$add_argument("-z", "--zscores", default = "3.5,5,10", metavar = "csv",
    help = "comma-separated list of Z-score thresholds [default %(default)s]")
parser$add_argument("-s", "--sexes", default = ",M,F", metavar = "csv",
    help = paste0("comma-separated list of sex filters; empty token = both ",
                  "[default '%(default)s']"))
parser$add_argument("-o", "--output_dir", default = "./", metavar = "directory",
    help = "directory containing regression CSV outputs [default %(default)s]")
parser$add_argument("--pval", type = "double", default = 0.05, metavar = "double",
    help = "significance threshold [default %(default)s]")
parser$add_argument("--fdr", action = "store_true",
    help = "apply BH FDR correction per contrast per condition")
parser$add_argument("--top_n", type = "integer", default = 20, metavar = "int",
    help = "number of top hits to display [default %(default)s]")
parser$add_argument("--skip_per_condition", action = "store_true",
    help = "skip rendering the 9 individual per-condition reports")
parser$add_argument("--per_condition_subdir", default = "per_condition_reports",
    metavar = "dir",
    help = paste0("subdirectory (relative to output_dir) for per-condition ",
                  "reports [default %(default)s]"))
parser$add_argument("--outfile", default = NULL, metavar = "filename",
    help = "output HTML filename (placed in output_dir) [default: auto-named]")
opt <- parser$parse_args()

if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("Package 'rmarkdown' is required. Install with: install.packages('rmarkdown')")

# ── locate Rmd files (must live alongside this script) ──────────────────────
combined_rmd     <- file.path(scriptdir, "Three_Group_Combined_Report.Rmd")
per_condition_rmd <- file.path(scriptdir, "Three_Group_Summary_Report.Rmd")
if (!file.exists(combined_rmd))
  stop(paste("Combined Rmd not found:", combined_rmd))
if (!opt$skip_per_condition && !file.exists(per_condition_rmd))
  stop(paste("Per-condition Rmd not found:", per_condition_rmd))

# ── parse condition lists ───────────────────────────────────────────────────
parse_csv <- function(s, allow_empty_tokens = FALSE) {
  parts <- trimws(strsplit(s, ",")[[1]])
  if (allow_empty_tokens) {
    if (grepl("^,", s)) parts <- c("", parts)
    if (grepl(",$", s)) parts <- c(parts, "")
    parts <- unique(parts)
  }
  parts
}
zscores <- parse_csv(opt$zscores)
sexes   <- parse_csv(opt$sexes, allow_empty_tokens = TRUE)

zscore_token <- function(z) {
  zn <- suppressWarnings(as.numeric(z))
  if (!is.na(zn) && zn == as.integer(zn)) as.character(as.integer(zn)) else z
}

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
per_cond_full <- file.path(opt$output_dir, opt$per_condition_subdir)

cat("Combined report — configuration:\n")
cat("  Groups   :", opt$group1, "/", opt$group2, "/", opt$group3, "\n")
cat("  Z-scores :", paste(zscores, collapse = ", "), "\n")
cat("  Sex      :", paste(ifelse(nchar(sexes) == 0, "(both)", sexes),
                          collapse = ", "), "\n")
cat("  Cond.    :", length(zscores) * length(sexes), "total\n")
cat("  Data dir :", opt$output_dir, "\n")
cat("  Per-cond.:", if (opt$skip_per_condition) "SKIPPED"
                   else per_cond_full, "\n")
cat("  FDR      :", opt$fdr, "\n\n")

# ── render per-condition reports (unless skipped) ───────────────────────────
if (!opt$skip_per_condition) {
  dir.create(per_cond_full, showWarnings = FALSE, recursive = TRUE)
  n_total <- length(zscores) * length(sexes)
  idx <- 0
  for (z in zscores) for (s in sexes) {
    idx <- idx + 1
    z_tok <- zscore_token(z)
    s_lbl <- if (nchar(s) == 0) "both" else s
    outfile <- paste0(
      gsub(" ", "_",
           paste("Three_Group_Summary_Report",
                 opt$group1, "vs", opt$group2, "vs", opt$group3,
                 "Z", z_tok, "sex", s_lbl, sep = "-")),
      ".html")
    outpath <- file.path(per_cond_full, outfile)

    cat(sprintf("[%d/%d] Per-condition: Z=%s sex=%s ... ",
                idx, n_total, z_tok, s_lbl))
    res <- tryCatch({
      rmarkdown::render(
        input       = per_condition_rmd,
        output_file = outpath,
        params      = list(
          group1      = opt$group1,
          group2      = opt$group2,
          group3      = opt$group3,
          zscore      = z_tok,
          sex         = s,
          output_dir  = opt$output_dir,
          pval_thresh = opt$pval,
          fdr         = opt$fdr,
          top_n       = opt$top_n
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )
      "OK"
    }, error = function(e) {
      paste("FAILED:", conditionMessage(e))
    })
    cat(res, "\n")
  }
  cat("\n")
}

# ── render combined report ──────────────────────────────────────────────────
if (is.null(opt$outfile)) {
  opt$outfile <- paste0(
    gsub(" ", "_",
         paste("Three_Group_Combined_Report",
               opt$group1, "vs", opt$group2, "vs", opt$group3,
               sep = "-")),
    ".html")
}
combined_out <- file.path(opt$output_dir, opt$outfile)

cat("Rendering combined report → ", combined_out, "\n")
rmarkdown::render(
  input       = combined_rmd,
  output_file = combined_out,
  params      = list(
    group1               = opt$group1,
    group2               = opt$group2,
    group3               = opt$group3,
    output_dir           = opt$output_dir,
    zscores              = paste(zscores, collapse = ","),
    sexes                = opt$sexes,
    pval_thresh          = opt$pval,
    fdr                  = opt$fdr,
    top_n                = opt$top_n,
    per_condition_subdir = opt$per_condition_subdir
  ),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

cat("\nDone.\n")
cat("Combined report:  ", combined_out, "\n")
if (!opt$skip_per_condition)
  cat("Per-condition dir:", per_cond_full, "\n")
