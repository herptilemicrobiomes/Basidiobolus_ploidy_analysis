#!/usr/bin/env Rscript
# summarize_coverage_quantized_lcg.R
#
# For each strain × reference combination in results/mosdepth_lcg/:
#   1. Load <strain>.<refbase>.quantized.bed.gz (BED4: chrom, start, end, label)
#   2. Produce per-strain multi-page PDFs (50 scaffolds per page):
#        heatmap - geom_rect tiles coloured by coverage class, all scaffolds stacked
#        bar     - stacked bar showing proportion of each class per scaffold
#
# Output: results/mosdepth_quantized_lcg/<refbase>/<strain>/{heatmap,bar}.pdf
#
# REFBASEs are auto-detected from files present in mosdepth_lcg/.
# Optionally restrict to one refbase by setting REFBASE env var:
#   Rscript summarize_coverage_quantized_lcg.R CBS931.73

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# ---- paths ----
samples_file <- "lcg.csv"
mosdepth_dir <- file.path("results", "mosdepth_lcg")
out_base     <- file.path("results", "mosdepth_quantized_lcg")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

SCAFFOLDS_PER_PAGE <- 50
MIN_SCAFFOLD_LENGTH <- 50000

# ---- load strain IDs from lcg.csv ----
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
strains <- samples[[1]]   # first column = SampleID

# ---- auto-detect available reference genomes ----
all_beds <- list.files(mosdepth_dir, pattern = "\\.quantized\\.bed\\.gz$")

# file stem = <strain>.<refbase>, strip the trailing .quantized.bed.gz
stems <- sub("\\.quantized\\.bed\\.gz$", "", all_beds)

# refbase = everything after the first component that matches a known strain
# e.g. "UHM1013.CBS931.73" -> strain="UHM1013", refbase="CBS931.73"
refbases_detected <- unique(sub(
  paste0("^(", paste(strains, collapse = "|"), ")\\."), "", stems
))
refbases_detected <- refbases_detected[refbases_detected != stems]  # drop non-matches

# allow CLI override: Rscript script.R CBS931.73
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  refbases <- args[1]
  message("Using REFBASE from command line: ", refbases)
} else {
  refbases <- refbases_detected
  message("Auto-detected REFBASEs: ", paste(refbases, collapse = ", "))
}

# ---- quantized coverage levels: cool -> hot ----
cov_levels <- c("NO_COVERAGE", "LOW_COVERAGE", "CALLABLE",
                "HIGH_COVERAGE", "VERY_HIGH_COVERAGE")

cov_colours <- c(
  NO_COVERAGE        = "#313695",
  LOW_COVERAGE       = "#74add1",
  CALLABLE           = "#ffffbf",
  HIGH_COVERAGE      = "#f46d43",
  VERY_HIGH_COVERAGE = "#a50026"
)

# ---- shared theme ----
base_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "#e8eaf6"),
    panel.spacing    = unit(0.3, "lines"),
    legend.position  = "right",
    axis.title       = element_text(size = 9)
  )

chunk <- function(x, n) split(x, ceiling(seq_along(x) / n))

# ===========================================================================
# Main loop: refbase × strain
# ===========================================================================
for (refbase in refbases) {
  message("\n====== REFBASE: ", refbase, " ======")
  ref_out <- file.path(out_base, refbase)
  dir.create(ref_out, showWarnings = FALSE, recursive = TRUE)

  for (strain in strains) {
    message("\n  === ", strain, " ===")

    bed_file <- file.path(mosdepth_dir,
                          paste0(strain, ".", refbase, ".quantized.bed.gz"))
    if (!file.exists(bed_file)) {
      message("  skipping (file not found): ", bed_file)
      next
    }

    bed <- read_tsv(bed_file,
                    col_names = c("chrom", "start", "end", "label"),
                    col_types = cols(
                      chrom = col_character(),
                      start = col_double(),
                      end   = col_double(),
                      label = col_character()
                    ))

    scaffold_order <- unique(bed$chrom)

    # filter short scaffolds
    scaffold_lengths <- bed |>
      group_by(chrom) |>
      summarise(length = max(end) - min(start), .groups = "drop")
    keep_scaffolds <- scaffold_lengths |>
      filter(length >= MIN_SCAFFOLD_LENGTH) |>
      pull(chrom)
    n_dropped      <- length(scaffold_order) - length(keep_scaffolds)
    scaffold_order <- scaffold_order[scaffold_order %in% keep_scaffolds]
    if (n_dropped > 0)
      message("  dropped ", n_dropped, " scaffolds shorter than ",
              MIN_SCAFFOLD_LENGTH, " bp")

    if (length(scaffold_order) == 0) {
      warning("  no scaffolds remain after length filter for ", strain)
      next
    }

    dat <- bed |>
      filter(chrom %in% scaffold_order) |>
      mutate(
        chrom = factor(chrom, levels = scaffold_order),
        label = factor(label, levels = cov_levels),
        width = end - start
      )

    nscaff      <- length(scaffold_order)
    out_dir     <- file.path(ref_out, strain)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    strain_label <- str_replace_all(strain, "\\.", " ")
    pages        <- chunk(scaffold_order, SCAFFOLDS_PER_PAGE)
    npages       <- length(pages)

    message("  ", nscaff, " scaffolds -> ", npages, " pages (ref: ", refbase, ")")

    # -----------------------------------------------------------------------
    # 1. HEATMAP
    # -----------------------------------------------------------------------
    heat_path <- file.path(out_dir, "heatmap.pdf")
    pdf(heat_path, width = 14, height = 20)

    for (i in seq_along(pages)) {
      scaffs   <- pages[[i]]
      page_dat <- dat |> filter(chrom %in% scaffs) |>
                    mutate(chrom = factor(chrom, levels = scaffs))

      p_heat <- ggplot(page_dat, aes(xmin = start, xmax = end,
                                     ymin = 0,    ymax = 1,
                                     fill = label)) +
        geom_rect() +
        scale_fill_manual(values = cov_colours, drop = FALSE,
                          name = "Coverage\nclass") +
        facet_wrap(~chrom, ncol = 1, scales = "free_x",
                   strip.position = "left") +
        labs(
          title    = paste0(strain_label, " vs ", refbase,
                            " - scaffold coverage (quantized heatmap)",
                            "  [page ", i, " / ", npages, "]"),
          subtitle = "Colour = mosdepth quantized coverage class",
          x        = "Genomic position (bp)",
          y        = NULL
        ) +
        base_theme +
        theme(
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          axis.text.y       = element_blank(),
          axis.ticks.y      = element_blank(),
          panel.grid        = element_blank()
        )

      print(p_heat)
    }

    dev.off()
    message("  saved -> ", heat_path)

    # -----------------------------------------------------------------------
    # 2. STACKED BAR
    # -----------------------------------------------------------------------
    bar_dat <- dat |>
      group_by(chrom, label) |>
      summarise(total_bp = sum(width), .groups = "drop") |>
      group_by(chrom) |>
      mutate(pct = total_bp / sum(total_bp) * 100) |>
      ungroup()

    bar_path <- file.path(out_dir, "bar.pdf")
    pdf(bar_path, width = 12, height = 18)

    for (i in seq_along(pages)) {
      scaffs       <- pages[[i]]
      page_bar_dat <- bar_dat |> filter(chrom %in% scaffs) |>
                        mutate(chrom = factor(chrom, levels = rev(scaffs)))

      p_bar <- ggplot(page_bar_dat, aes(x = pct, y = chrom, fill = label)) +
        geom_col(width = 0.8) +
        scale_fill_manual(values = cov_colours, drop = FALSE,
                          name = "Coverage\nclass") +
        scale_x_continuous(expand = c(0, 0),
                           labels = function(x) paste0(x, "%")) +
        labs(
          title    = paste0(strain_label, " vs ", refbase,
                            " - scaffold coverage class proportions",
                            "  [page ", i, " / ", npages, "]"),
          subtitle = "Percentage of scaffold bases in each quantized class",
          x        = "Percentage of scaffold (bp)",
          y        = "Scaffold"
        ) +
        base_theme

      print(p_bar)
    }

    dev.off()
    message("  saved -> ", bar_path)
  }
}

message("\nAll done. Output in: ", out_base)
