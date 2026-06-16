#!/usr/bin/env Rscript
# Summarize mosdepth coverage data across strains
# Reads all results/mosdepth/*.mosdepth.summary.txt files
# Produces histogram and boxplot of per-scaffold coverage depth,
# distinguishing _region scaffolds from whole-scaffold entries.

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ---- locate files ----
mosdepth_dir <- file.path("results", "mosdepth_lcg")
summary_files <- list.files(mosdepth_dir,
                            pattern = "\\.mosdepth\\.summary\\.txt$",
                            full.names = TRUE)

if (length(summary_files) == 0) {
  stop("No *.mosdepth.summary.txt files found in ", mosdepth_dir)
}

# ---- read and combine ----
read_summary <- function(path) {
  strain <- str_remove(basename(path), "\\.mosdepth\\.summary\\.txt$")
  df <- read_tsv(path, col_types = cols(
    chrom  = col_character(),
    length = col_double(),
    bases  = col_double(),
    mean   = col_double(),
    min    = col_double(),
    max    = col_double()
  ))
  df$strain <- strain
  df
}

dat <- bind_rows(lapply(summary_files, read_summary))

# ---- filter and classify ----
dat <- dat |>
  filter(!chrom %in% c("total", "total_region")) |>
  mutate(type = if_else(str_ends(chrom, "_region"), "region", "whole_scaffold"))

# column 4 (mean) is the depth of coverage
# dat$mean is already that column

# ---- output directory ----
out_dir <- file.path("results", "mosdepth_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- colour palette ----
type_colors <- c(whole_scaffold = "#2166ac", region = "#d6604d")
type_labels  <- c(whole_scaffold = "Whole scaffold", region = "Region (BED)")

# ---- median-normalise within each strain (for ploidy comparison) ----
# Values > 1 indicate higher-than-median copy number; 2.0 = doubled coverage
dat <- dat |>
  group_by(strain) |>
  mutate(median_cov   = median(mean, na.rm = TRUE),
         mean_norm    = mean / median_cov) |>
  ungroup()

# ===========================================================================
# 1. Histogram – log10 x-axis, one panel per strain
#    log10 scale: each unit = one doubling; ploidy levels appear evenly spaced
# ===========================================================================
p_hist <- ggplot(dat, aes(x = mean, fill = type)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.65, colour = NA) +
  facet_wrap(~strain, scales = "free_y") +
  scale_x_continuous(
    trans  = "log10",
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_fill_manual(values = type_colors, labels = type_labels,
                    name = "Scaffold type") +
  labs(
    title = "Distribution of per-scaffold mean coverage depth by strain (log10 scale)",
    x     = "Mean depth of coverage (log10)",
    y     = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

hist_out <- file.path(out_dir, "coverage_histogram_LCG.pdf")
ggsave(hist_out, p_hist, width = 14, height = 10)
message("Saved histogram -> ", hist_out)

# ===========================================================================
# 2. Box-and-whisker – log10 y-axis, one box per scaffold type, panels per strain
# ===========================================================================
p_box <- ggplot(dat, aes(x = type, y = mean, fill = type)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.5) +
  facet_wrap(~strain, scales = "free_y") +
  scale_y_continuous(
    trans  = "log10",
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_fill_manual(values = type_colors, labels = type_labels,
                    name = "Scaffold type") +
  scale_x_discrete(labels = type_labels) +
  labs(
    title = "Per-scaffold mean coverage depth by strain and scaffold type (log10 scale)",
    x     = "Scaffold type",
    y     = "Mean depth of coverage (log10)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    legend.position  = "none"
  )

box_out <- file.path(out_dir, "coverage_boxplot_LCG.pdf")
ggsave(box_out, p_box, width = 14, height = 10)
message("Saved boxplot  -> ", box_out)

# ===========================================================================
# 3. Median-normalised boxplot – all strains on one shared y-axis
#    Each strain is normalised by its own median coverage so that 1.0 = median,
#    2.0 = doubled coverage (potential diploid vs haploid shift), etc.
#    This removes sequencing-depth differences and makes strains directly comparable.
# ===========================================================================
p_norm <- ggplot(dat, aes(x = strain, y = mean_norm, fill = type)) +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6,
               position = position_dodge(width = 0.7)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(0.25, 0.5, 1, 1.5, 2, 3, 4),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_fill_manual(values = type_colors, labels = type_labels,
                    name = "Scaffold type") +
  labs(
    title = "Median-normalised coverage across strains (log10 scale)",
    subtitle = "1.0 = strain median; dashed lines at 0.5×, 1×, 1.5×, 2× median",
    x     = "Strain",
    y     = "Coverage / strain median (log10)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "bottom"
  )

norm_out <- file.path(out_dir, "coverage_normalised_boxplot_LCG.pdf")
ggsave(norm_out, p_norm, width = 14, height = 8)
message("Saved normalised boxplot -> ", norm_out)
