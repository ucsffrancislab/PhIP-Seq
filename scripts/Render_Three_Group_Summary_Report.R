#!/usr/bin/env Rscript

#   Render_Three_Group_Summary_Report.R
#
#   Renders Three_Group_Summary_Report.Rmd into a self-contained HTML file.
#   The Rmd must live in the same directory as this script (alongside GenoLib.R
#   and the three Multi_Plate_Three_Group_*.R scripts).
#
#   Usage:
#       module load r
#       Render_Three_Group_Summary_Report.R \
#           -a "PF Patient" \
#           -b "Endemic Control" \
#           -c "Non Endemic Control" \
#           -o /path/to/output_dir \
#           [--pval 0.05] [--fdr] [--top_n 20] [--outfile report.html]
#
#   The script searches output_dir for CSV files matching the naming patterns
#   produced by:
#       Multi_Plate_Three_Group_Peptide_Regression.R
#       Multi_Plate_Three_Group_VirScan_Seropositivity_Regression.R
#       Multi_Plate_Three_Group_VirHitFrac_Seropositivity_Regression.R
#   It will use the most recently modified matching file for each type.
#
#   Dependencies: rmarkdown, knitr (and all packages listed in the Rmd)

library("argparse")

args       <- commandArgs()
scriptdir  <- dirname(sub("--file=", "", args[grepl("--file=", args)]))

parser <- ArgumentParser()
parser$add_argument("-a", "--group1",  required = TRUE, metavar = "group",
    help = "first group label  (reference / cases)")
parser$add_argument("-b", "--group2",  required = TRUE, metavar = "group",
    help = "second group label (control 1)")
parser$add_argument("-c", "--group3",  required = TRUE, metavar = "group",
    help = "third group label  (control 2)")
parser$add_argument("-o", "--output_dir", default = "./", metavar = "directory",
    help = "directory containing regression CSV outputs [default %(default)s]")
parser$add_argument("--pval", type = "double", default = 0.05, metavar = "double",
    help = "significance threshold [default %(default)s]")
parser$add_argument("--fdr", action = "store_true",
    help = "apply BH FDR correction before thresholding")
parser$add_argument("--top_n", type = "integer", default = 20, metavar = "integer",
    help = "number of top hits to display per section [default %(default)s]")
parser$add_argument("--outfile", default = NULL, metavar = "filename",
    help = "output HTML filename (placed in output_dir) [default: auto-named]")
opt <- parser$parse_args()

if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("Package 'rmarkdown' is required. Install with: install.packages('rmarkdown')")

rmd_path <- file.path(scriptdir, "Three_Group_Summary_Report.Rmd")
if (!file.exists(rmd_path))
  stop(paste("Rmd not found:", rmd_path))

# Auto-name the output file if not specified
if (is.null(opt$outfile)) {
  opt$outfile <- paste0(
    gsub(" ", "_",
         paste("Three_Group_Summary_Report",
               opt$group1, "vs", opt$group2, "vs", opt$group3,
               sep = "-")),
    ".html")
}

outfile_path <- file.path(opt$output_dir, opt$outfile)
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Rendering report...\n")
cat("  Groups  :", opt$group1, "/", opt$group2, "/", opt$group3, "\n")
cat("  Data dir:", opt$output_dir, "\n")
cat("  Output  :", outfile_path, "\n")
cat("  FDR     :", opt$fdr, "\n")
cat("  p thresh:", opt$pval, "\n\n")

rmarkdown::render(
  input       = rmd_path,
  output_file = outfile_path,
  params      = list(
    group1       = opt$group1,
    group2       = opt$group2,
    group3       = opt$group3,
    output_dir   = opt$output_dir,
    pval_thresh  = opt$pval,
    fdr          = opt$fdr,
    top_n        = opt$top_n
  ),
  envir  = new.env(parent = globalenv()),
  quiet  = FALSE
)

cat("\nReport written to:", outfile_path, "\n")
