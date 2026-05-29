#!/usr/bin/env Rscript

#   Three_Group_Volcano_Grid.R
#
#   Produces a 2×2 grid of volcano-style plots from the output of any of the
#   three Multi_Plate_Three_Group_*.R regression scripts:
#
#       Multi_Plate_Three_Group_Peptide_Regression.R
#       Multi_Plate_Three_Group_VirHitFrac_Seropositivity_Regression.R
#       Multi_Plate_Three_Group_VirScan_Seropositivity_Regression.R
#
#   The four panels are:
#
#       [1] Side-by-side dual volcanos (two panels drawn as one, g2-vs-g1 left,
#           g3-vs-g1 right) shown as a single combined top row.
#           Actually implemented as two separate ggplot objects assembled with
#           patchwork, sharing a common y-axis range.
#
#       Top-left  (panel 1): Volcano — group2 vs. group1
#                            x = beta_<g2>_vs_<g1>, y = -log10(pval_<g2>_vs_<g1>)
#
#       Top-right (panel 2): Volcano — group3 vs. group1
#                            x = beta_<g3>_vs_<g1>, y = -log10(pval_<g3>_vs_<g1>)
#
#       Bottom-left  (panel 3): Mirrored / back-to-back volcano
#                            Both contrasts share the x-axis (beta).
#                            g2-vs-g1 points plotted upward (+log10 p),
#                            g3-vs-g1 points plotted downward (-log10 p).
#
#       Bottom-right (panel 4): Contrast scatter
#                            x = beta_<g2>_vs_<g1>,
#                            y = beta_<g3>_vs_<g1>,
#                            colour = -log10(min p-value across contrasts),
#                            size   = -log10(min p-value).
#
#   Points are coloured by significance category in panels 1–3:
#       "Both"    : significant in both contrasts
#       "G2 only" : significant in g2-vs-g1 contrast only
#       "G3 only" : significant in g3-vs-g1 contrast only
#       "Neither" : not significant in either
#
#   Usage:
#       module load r
#       Three_Group_Volcano_Grid.R \
#           --input /path/to/Multiplate_Three_Group_Peptide_Comparison_...csv \
#           --group1 "PF Patient" \
#           --group2 "Non Endemic Control" \
#           --group3 "Healthy Control" \
#           --label_col species \
#           --pval_threshold 0.05 \
#           --output_dir /path/to/output \
#           --output_basename my_volcano_grid
#
#   Optional flags:
#       --fdr              apply BH FDR correction before thresholding
#       --label_top N      label the top N significant points (default 10)
#       --label_col COL    column to use for point labels (default: species;
#                          use 'peptide' for peptide-level output)
#       --width  W         plot width in inches  (default 14)
#       --height H         plot height in inches (default 12)
#
#   Dependencies: ggplot2, patchwork, ggrepel, scales, data.table, argparse

options(warn = 1)

# ---------------------------------------------------------------------------
# Package loading with informative error
# ---------------------------------------------------------------------------

load_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste0("Package '", pkg, "' is required but not installed.\n",
                    "Install with: install.packages('", pkg, "')"), call. = FALSE)
    }
    library(pkg, character.only = TRUE)
}

load_pkg("argparse")
load_pkg("data.table")
load_pkg("ggplot2")
load_pkg("patchwork")
load_pkg("ggrepel")
load_pkg("scales")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

args       <- commandArgs()
scriptname <- sub("--file=", "", args[grepl("--file=", args)])

parser <- ArgumentParser(description = scriptname)
parser$add_argument("-i", "--input", type = "character", required = TRUE,
    help = "path to the three-group regression output CSV",
    metavar = "CSV file")
parser$add_argument("-a", "--group1", type = "character", required = TRUE,
    help = "reference group label (must match column names in CSV)",
    metavar = "group")
parser$add_argument("-b", "--group2", type = "character", required = TRUE,
    help = "second group label", metavar = "group")
parser$add_argument("-c", "--group3", type = "character", required = TRUE,
    help = "third group label", metavar = "group")
parser$add_argument("--label_col", type = "character", default = "species",
    help = "column to use for point labels [default %(default)s]",
    metavar = "column")
parser$add_argument("--pval_threshold", type = "double", default = 0.05,
    help = "significance threshold (pre- or post-FDR) [default %(default)s]",
    metavar = "double")
parser$add_argument("--fdr", action = "store_true",
    help = "apply Benjamini-Hochberg FDR correction before thresholding")
parser$add_argument("--label_top", type = "integer", default = 10,
    help = "number of top significant points to label [default %(default)s]",
    metavar = "integer")
parser$add_argument("-o", "--output_dir", type = "character", default = "./",
    help = "output directory [default %(default)s]", metavar = "directory")
parser$add_argument("--output_basename", type = "character",
    default = "three_group_volcano_grid",
    help = "output file basename without extension [default %(default)s]",
    metavar = "basename")
parser$add_argument("--width",  type = "double", default = 14,
    help = "plot width in inches  [default %(default)s]", metavar = "inches")
parser$add_argument("--height", type = "double", default = 12,
    help = "plot height in inches [default %(default)s]", metavar = "inches")
opt <- parser$parse_args()

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Derive column names from group labels (mirrors regression script logic)
# ---------------------------------------------------------------------------

g2_vs_g1 <- paste0(opt$group2, "_vs_", opt$group1)
g3_vs_g1 <- paste0(opt$group3, "_vs_", opt$group1)

col_beta_g2  <- paste0("beta_", g2_vs_g1)
col_se_g2    <- paste0("se_",   g2_vs_g1)
col_pval_g2  <- paste0("pval_", g2_vs_g1)
col_beta_g3  <- paste0("beta_", g3_vs_g1)
col_se_g3    <- paste0("se_",   g3_vs_g1)
col_pval_g3  <- paste0("pval_", g3_vs_g1)

# ---------------------------------------------------------------------------
# Read and validate input
# ---------------------------------------------------------------------------

print(paste("Reading:", opt$input))
dat <- data.frame(data.table::fread(opt$input, sep = ",", header = TRUE),
                  check.names = FALSE)

required_cols <- c(col_beta_g2, col_pval_g2, col_beta_g3, col_pval_g3)
missing <- setdiff(required_cols, colnames(dat))
if (length(missing) > 0) {
    stop(paste("Missing expected columns:\n",
               paste(missing, collapse = "\n"),
               "\nCheck that --group1/2/3 match the labels used when running the regression."),
         call. = FALSE)
}

if (!(opt$label_col %in% colnames(dat))) {
    warning(paste0("--label_col '", opt$label_col,
                   "' not found; falling back to row numbers."))
    dat$label_col_used <- as.character(seq_len(nrow(dat)))
    opt$label_col <- "label_col_used"
}

# ---------------------------------------------------------------------------
# Coerce to numeric, drop rows where both p-values are NA
# ---------------------------------------------------------------------------

dat[[col_beta_g2]] <- as.numeric(dat[[col_beta_g2]])
dat[[col_pval_g2]] <- as.numeric(dat[[col_pval_g2]])
dat[[col_beta_g3]] <- as.numeric(dat[[col_beta_g3]])
dat[[col_pval_g3]] <- as.numeric(dat[[col_pval_g3]])

n_raw <- nrow(dat)
dat   <- dat[!(is.na(dat[[col_pval_g2]]) & is.na(dat[[col_pval_g3]])), ]
print(paste(n_raw - nrow(dat), "rows dropped (both p-values NA);",
            nrow(dat), "rows retained."))

# ---------------------------------------------------------------------------
# Optional FDR correction (applied per contrast independently)
# ---------------------------------------------------------------------------

if (opt$fdr) {
    print("Applying BH FDR correction per contrast.")
    dat[[col_pval_g2]] <- p.adjust(dat[[col_pval_g2]], method = "BH")
    dat[[col_pval_g3]] <- p.adjust(dat[[col_pval_g3]], method = "BH")
}

# ---------------------------------------------------------------------------
# Derived columns used across all four panels
# ---------------------------------------------------------------------------

alpha <- opt$pval_threshold

dat$neglog10_p_g2  <- -log10(dat[[col_pval_g2]])
dat$neglog10_p_g3  <- -log10(dat[[col_pval_g3]])
dat$pval_min       <- pmin(dat[[col_pval_g2]], dat[[col_pval_g3]], na.rm = TRUE)
dat$neglog10_p_min <- -log10(dat$pval_min)

sig_g2 <- !is.na(dat[[col_pval_g2]]) & dat[[col_pval_g2]] < alpha
sig_g3 <- !is.na(dat[[col_pval_g3]]) & dat[[col_pval_g3]] < alpha

dat$sig_category <- ifelse(sig_g2 & sig_g3, "Both",
                    ifelse(sig_g2,            paste0(opt$group2, " only"),
                    ifelse(sig_g3,            paste0(opt$group3, " only"),
                                              "Neither")))
dat$sig_category <- factor(dat$sig_category,
    levels = c("Both",
               paste0(opt$group2, " only"),
               paste0(opt$group3, " only"),
               "Neither"))

# Significance colour palette — consistent across panels 1–3
sig_colours <- c(
    "Both"                          = "#C0392B",   # red
    "Neither"                       = "#BDC3C7"    # light grey
)
sig_colours[paste0(opt$group2, " only")] <- "#2980B9"   # blue
sig_colours[paste0(opt$group3, " only")] <- "#27AE60"   # green

# Points to label: top N by minimum p-value, among those significant in either
sig_any <- sig_g2 | sig_g3
if (sum(sig_any) > 0) {
    top_idx <- order(dat$pval_min[sig_any])[seq_len(min(opt$label_top, sum(sig_any)))]
    label_ids <- which(sig_any)[top_idx]
} else {
    label_ids <- integer(0)
}
dat$label_text <- ifelse(seq_len(nrow(dat)) %in% label_ids,
                          as.character(dat[[opt$label_col]]), NA_character_)

# Shared horizontal threshold line height
sig_line <- -log10(alpha)

# ---------------------------------------------------------------------------
# Theme shared across all panels
# ---------------------------------------------------------------------------

theme_volcano <- function(base_size = 11) {
    theme_bw(base_size = base_size) +
    theme(
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = "grey92"),
        plot.title        = element_text(face = "bold", size = base_size + 1,
                                         hjust = 0.5),
        plot.subtitle     = element_text(hjust = 0.5, colour = "grey40",
                                         size = base_size - 1),
        legend.position   = "bottom",
        legend.title      = element_text(face = "bold"),
        legend.key.size   = unit(0.4, "cm"),
        axis.title        = element_text(size = base_size)
    )
}

# ---------------------------------------------------------------------------
# Panel 1 — Volcano: group2 vs. group1
# ---------------------------------------------------------------------------

xlim1 <- range(dat[[col_beta_g2]], na.rm = TRUE)
xlim1 <- c(-max(abs(xlim1)), max(abs(xlim1)))  # symmetric

p1 <- ggplot(dat, aes(x = .data[[col_beta_g2]],
                      y = neglog10_p_g2,
                      colour = sig_category,
                      label = label_text)) +
    geom_hline(yintercept = sig_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid",
               colour = "grey70", linewidth = 0.4) +
    geom_point(alpha = 0.65, size = 1.6) +
    ggrepel::geom_text_repel(
        na.rm = TRUE, size = 2.8, max.overlaps = 20,
        segment.colour = "grey50", segment.size = 0.3,
        show.legend = FALSE) +
    scale_colour_manual(values = sig_colours, name = "Significant in",
                        drop = FALSE) +
    scale_x_continuous(limits = xlim1, oob = scales::squish) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
        title    = paste0(opt$group2, " vs. ", opt$group1),
        subtitle = if (opt$fdr) paste0("BH FDR < ", alpha) else paste0("p < ", alpha),
        x        = paste0("Beta (", opt$group2, " vs. ", opt$group1, ")"),
        y        = expression(-log[10](p))
    ) +
    theme_volcano()

# ---------------------------------------------------------------------------
# Panel 2 — Volcano: group3 vs. group1
# ---------------------------------------------------------------------------

xlim2 <- range(dat[[col_beta_g3]], na.rm = TRUE)
xlim2 <- c(-max(abs(xlim2)), max(abs(xlim2)))

p2 <- ggplot(dat, aes(x = .data[[col_beta_g3]],
                      y = neglog10_p_g3,
                      colour = sig_category,
                      label = label_text)) +
    geom_hline(yintercept = sig_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid",
               colour = "grey70", linewidth = 0.4) +
    geom_point(alpha = 0.65, size = 1.6) +
    ggrepel::geom_text_repel(
        na.rm = TRUE, size = 2.8, max.overlaps = 20,
        segment.colour = "grey50", segment.size = 0.3,
        show.legend = FALSE) +
    scale_colour_manual(values = sig_colours, name = "Significant in",
                        drop = FALSE) +
    scale_x_continuous(limits = xlim2, oob = scales::squish) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
        title    = paste0(opt$group3, " vs. ", opt$group1),
        subtitle = if (opt$fdr) paste0("BH FDR < ", alpha) else paste0("p < ", alpha),
        x        = paste0("Beta (", opt$group3, " vs. ", opt$group1, ")"),
        y        = expression(-log[10](p))
    ) +
    theme_volcano()

# ---------------------------------------------------------------------------
# Panel 3 — Mirrored / back-to-back volcano
#
# g2-vs-g1 plotted upward   (positive  -log10 p)
# g3-vs-g1 plotted downward (negative  -log10 p)
# Both share the same x-axis (beta) with a symmetric range.
# ---------------------------------------------------------------------------

mirror_dat <- rbind(
    data.frame(
        beta      = dat[[col_beta_g2]],
        mirror_y  = dat$neglog10_p_g2,    # positive → upward
        sig_cat   = dat$sig_category,
        label_txt = dat$label_text,
        contrast  = paste0(opt$group2, " vs. ", opt$group1)
    ),
    data.frame(
        beta      = dat[[col_beta_g3]],
        mirror_y  = -dat$neglog10_p_g3,   # negative → downward
        sig_cat   = dat$sig_category,
        label_txt = dat$label_text,
        contrast  = paste0(opt$group3, " vs. ", opt$group1)
    )
)
mirror_dat <- mirror_dat[!is.na(mirror_dat$beta), ]

xlim3    <- range(mirror_dat$beta, na.rm = TRUE)
xlim3    <- c(-max(abs(xlim3)), max(abs(xlim3)))
ylim3    <- range(mirror_dat$mirror_y, na.rm = TRUE)
ylim3abs <- max(abs(ylim3))

# Custom y-axis: show absolute value but label direction
y_breaks <- pretty(c(-ylim3abs, ylim3abs), n = 8)
y_labels <- as.character(abs(y_breaks))

p3 <- ggplot(mirror_dat, aes(x = beta, y = mirror_y,
                              colour = sig_cat, label = label_txt)) +
    # centre mirror line
    geom_hline(yintercept = 0,         colour = "grey30", linewidth = 0.6) +
    # significance threshold lines (both directions)
    geom_hline(yintercept =  sig_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -sig_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0,         colour = "grey70", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.4) +
    ggrepel::geom_text_repel(
        data = mirror_dat[!is.na(mirror_dat$label_txt), ],
        na.rm = TRUE, size = 2.5, max.overlaps = 20,
        segment.colour = "grey50", segment.size = 0.3,
        show.legend = FALSE) +
    # group direction annotations
    annotate("text", x = xlim3[1] * 0.95, y =  ylim3abs * 0.92,
             label = paste0("\u2191 ", opt$group2, " vs. ", opt$group1),
             hjust = 0, size = 3, colour = "grey30", fontface = "italic") +
    annotate("text", x = xlim3[1] * 0.95, y = -ylim3abs * 0.92,
             label = paste0("\u2193 ", opt$group3, " vs. ", opt$group1),
             hjust = 0, size = 3, colour = "grey30", fontface = "italic") +
    scale_colour_manual(values = sig_colours, name = "Significant in",
                        drop = FALSE) +
    scale_x_continuous(limits = xlim3, oob = scales::squish) +
    scale_y_continuous(limits = c(-ylim3abs, ylim3abs),
                       breaks = y_breaks, labels = y_labels,
                       expand = expansion(mult = 0.05)) +
    labs(
        title    = "Mirrored volcano",
        subtitle = paste0(opt$group2, " \u2191  |  ", opt$group3, " \u2193"),
        x        = "Beta",
        y        = expression(-log[10](p))
    ) +
    theme_volcano()

# ---------------------------------------------------------------------------
# Panel 4 — Contrast scatter
#
# x = beta_g2_vs_g1,  y = beta_g3_vs_g1
# colour = -log10(min p-value), size = -log10(min p-value)
# Quadrant lines at x=0 and y=0; significance ellipses optional.
# ---------------------------------------------------------------------------

scatter_dat <- dat[!is.na(dat[[col_beta_g2]]) & !is.na(dat[[col_beta_g3]]), ]

xlim4 <- range(scatter_dat[[col_beta_g2]], na.rm = TRUE)
xlim4 <- c(-max(abs(xlim4)), max(abs(xlim4)))
ylim4 <- range(scatter_dat[[col_beta_g3]], na.rm = TRUE)
ylim4 <- c(-max(abs(ylim4)), max(abs(ylim4)))

# Cap neglog10 for colour scale so a handful of extreme points
# don't compress all the others into the bottom of the palette
neglog10_cap <- quantile(scatter_dat$neglog10_p_min,
                         probs = 0.995, na.rm = TRUE)
scatter_dat$neglog10_p_min_capped <- pmin(scatter_dat$neglog10_p_min,
                                          neglog10_cap)

p4 <- ggplot(scatter_dat,
             aes(x     = .data[[col_beta_g2]],
                 y     = .data[[col_beta_g3]],
                 colour = neglog10_p_min_capped,
                 size   = neglog10_p_min_capped,
                 label  = label_text)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
    # diagonal (y=x) reference — features behaving identically in both contrasts
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey55", linewidth = 0.5) +
    geom_point(alpha = 0.70) +
    ggrepel::geom_text_repel(
        na.rm = TRUE, size = 2.8, max.overlaps = 20,
        segment.colour = "grey50", segment.size = 0.3,
        show.legend = FALSE) +
    scale_colour_gradientn(
        colours  = c("grey85", "#F39C12", "#E74C3C", "#8E44AD"),
        name     = expression(-log[10](p[min])),
        limits   = c(0, neglog10_cap),
        oob      = scales::squish
    ) +
    scale_size_continuous(
        range  = c(0.8, 4),
        limits = c(0, neglog10_cap),
        guide  = "none"          # size redundant with colour
    ) +
    scale_x_continuous(limits = xlim4, oob = scales::squish) +
    scale_y_continuous(limits = ylim4, oob = scales::squish) +
    # quadrant annotations
    annotate("text", x = xlim4[2] * 0.97, y = ylim4[2] * 0.97,
             label = "Both\nenriched", hjust = 1, vjust = 1,
             size = 2.8, colour = "grey40", fontface = "italic") +
    annotate("text", x = xlim4[1] * 0.97, y = ylim4[1] * 0.97,
             label = "Both\ndepleted", hjust = 0, vjust = 0,
             size = 2.8, colour = "grey40", fontface = "italic") +
    annotate("text", x = xlim4[2] * 0.97, y = ylim4[1] * 0.97,
             label = paste0(opt$group2, "\nenriched only"),
             hjust = 1, vjust = 0,
             size = 2.8, colour = "#2980B9", fontface = "italic") +
    annotate("text", x = xlim4[1] * 0.97, y = ylim4[2] * 0.97,
             label = paste0(opt$group3, "\nenriched only"),
             hjust = 0, vjust = 1,
             size = 2.8, colour = "#27AE60", fontface = "italic") +
    labs(
        title    = "Contrast scatter",
        subtitle = paste0("Diagonal = concordant effect in both contrasts"),
        x        = paste0("Beta (", opt$group2, " vs. ", opt$group1, ")"),
        y        = paste0("Beta (", opt$group3, " vs. ", opt$group1, ")")
    ) +
    theme_volcano() +
    theme(legend.position = "right")

# ---------------------------------------------------------------------------
# Assemble 2×2 grid with patchwork
# ---------------------------------------------------------------------------

# Collect the legend from panels 1–3 (they share a categorical colour scale)
# and place it below the grid.  Panel 4 keeps its own gradient legend on
# the right, handled by patchwork's guides collection.

grid <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
        title   = paste0("Three-group comparison: ", opt$group1,
                          " vs. ", opt$group2, " vs. ", opt$group3),
        subtitle = basename(opt$input),
        theme   = theme(
            plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40")
        )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

outfile_pdf <- file.path(opt$output_dir,
                         paste0(opt$output_basename, ".pdf"))
outfile_png <- file.path(opt$output_dir,
                         paste0(opt$output_basename, ".png"))

print(paste("Writing PDF:", outfile_pdf))
ggsave(outfile_pdf, plot = grid,
       width = opt$width, height = opt$height, device = "pdf")

print(paste("Writing PNG:", outfile_png))
ggsave(outfile_png, plot = grid,
       width = opt$width, height = opt$height, device = "png", dpi = 300)

print("Done.")
